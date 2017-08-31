/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.ngs;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.CovariateEnum;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.reader.SequencesReader;
import com.rtg.sam.BamIndexMerge;
import com.rtg.sam.BamIndexer;
import com.rtg.tabix.IndexingStreamCreator;
import com.rtg.tabix.TabixIndexMerge;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.IORunnable;
import com.rtg.util.IORunnableProxy;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.OneShotTimer;
import com.rtg.util.intervals.Range;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.util.io.FileUtils;

import htsjdk.samtools.SAMFileHeader;

/**
 * Multicore version of filter concatenation. This class determines
 * whether input and output files are compressed and uses appropriate
 * strategies for dividing the work among multiple threads and
 * concatenating.
 *
 */
@TestClass(value = {"com.rtg.ngs.TopNPairedEndOutputProcessorSyncTest"})
public abstract class AbstractMulticoreFilterConcat {

  // Each thread can zip 30MB/sec to end up writing 5MB/sec to disk. Disks can write ~50MB, so at most ~10 threads
  protected static final int MAX_FILTERCONCAT_THREADS = 10; //Integer.parseInt(System.getProperty("rtg.max_filterconcat_threads", "10"));

  protected final NgsParams mParams;

  protected String mThreadNamePrefix = "";

  /**
   * Creates a new <code>AbstractMulticoreFilterConcat</code> instance.
   *
   * @param params a <code>NgsParams</code> value containing
   * configuration information (such as the number of threads and
   * whether input and output files are to be compressed).
   */
  public AbstractMulticoreFilterConcat(NgsParams params) {
    mParams = params;
  }

  /**
   * Creates a new <code>AbstractMulticoreFilterConcat</code> instance.
   *
   * @param params a <code>NgsParams</code> value containing
   * configuration information (such as the number of threads and
   * whether input and output files are to be compressed).
   * @param threadNamePrefix prefix for thread name during logging (to differentiate different phases' filter concatenation jobs)
   */
  public AbstractMulticoreFilterConcat(NgsParams params, String threadNamePrefix) {
    mParams = params;
    mThreadNamePrefix = threadNamePrefix;
  }

  /**
   * Perform any post-processing on the intermediate files (read, re-write) that is required.
   * @param numThreads the number of threads (same as the number of intermediate files).
   * @param intermediate the intermediate files to be post-processed.
   * @param intermediateIndexes the array to populate with the new index files produced in any post-processing.
   * @param header the SAM file header.
   * @param createIndex whether indexes are being created
   * @throws IOException if an error occurs.
   */
  protected void postProcessIntermediateFiles(final int numThreads, final File[] intermediate, final File[] intermediateIndexes, SAMFileHeader header, boolean createIndex) throws IOException {
    // Default implementation does nothing
  }

  /**
   * Create and return output streams and index runners for given intermediate file.
   * @param numThreads number of threads being used to output intermediate files.
   * @param intermediate array of the intermediate files being created.
   * @param intermediateIndexes array of intermediate index files (must be populated as part of this method).
   * @param samGzipIntFiles if output SAM files should be compressed.
   * @param createIndex true if indexes should be created
   * @param i the number of the intermediate file to create the streams for.
   * @return an output wrapper containing the output stream for the intermediate file and the index runner for the index.
   * @throws IOException if an error occurs.
   */
  protected OutputWrapper createStreams(final int numThreads, final File[] intermediate, final File[] intermediateIndexes, final boolean samGzipIntFiles, boolean createIndex, int i) throws IOException {
    final OutputWrapper outWrapper;
    final File dataFile = intermediate[i];
    if (createIndex && SamSingleEndOutputProcessor.canIndex(mParams, intermediate[0])) {
      final PipedInputStream pipeToIndexIn = new PipedInputStream(); //closed by IndexRunner
      final PipedOutputStream pipeToIndexOut = new PipedOutputStream(pipeToIndexIn); //closed by SubFilter via intStream
      final OutputStream intStream = FileUtils.createTeedOutputStream(dataFile, pipeToIndexOut, samGzipIntFiles, false, samGzipIntFiles && (i == numThreads - 1));
      final File indexFilename = indexFileName(dataFile, mParams.outputParams().bam());
      intermediateIndexes[i] = indexFilename;
      final FileOutputStream indexOut = new FileOutputStream(indexFilename);
      final TabixIndexer.IndexerFactory indexer = mParams.outputParams().bam() ? null : new TabixIndexer.SamIndexerFactory();
      final IndexingStreamCreator.IndexRunner indexRunner = new IndexingStreamCreator.IndexRunner(pipeToIndexIn, indexOut, indexer, i == 0, (int) mParams.searchParams().numberSequences(), dataFile.toString());
      outWrapper = new OutputWrapper(intStream, indexRunner);
    } else {
      outWrapper = new OutputWrapper(FileUtils.createOutputStream(dataFile, samGzipIntFiles, samGzipIntFiles && (i == numThreads - 1)), null);
    }
    return outWrapper;
  }

  static File indexFileName(File dataFile, boolean bam) {
    if (bam) {
      return BamIndexer.indexFileName(dataFile);
    } else {
      return TabixIndexer.indexFileName(dataFile);
    }
  }

  /**
   * Filters and concatenates the supplied temp files to the supplied
   * out file using multiple threads.
   *
   * @param tempFiles the <code>Files</code> to be filtered and concatenated
   * @param outFile the destination <code>File</code>
   * @param header contains SAM header structure
   * @param outputParams parameters describing calibration settings
   * @return list of files at current stage of processing. Should be single SAM/BAM file unless <code>delayMerge</code> is set, in which case it is the filtered SAM/BAM files
   * @exception IOException if an error occurs.
   */
  FilterConcatIntermediateFiles filterConcat(File[] tempFiles, File outFile, SAMFileHeader header, NgsOutputParams outputParams) throws IOException {
    // @param delayMerge set this to prevent indexing and final merging of filtered files
    final boolean delayMerge = outputParams.unify();
    final boolean calibrate = outputParams.calibrate();
    final boolean samGzipIntFiles = mParams.outputParams().isCompressOutput();
    final boolean noLongSequences;
    if (mParams.searchParams().reader().maxLength() > TabixIndexer.MAXIMUM_REFERENCE_LENGTH) {
      Diagnostic.warning("Cannot produce TABIX index as maximum reference sequence length is exceeded.");
      noLongSequences = false;
    } else {
      noLongSequences = true;
    }
    final boolean createIndex = mParams.outputParams().outputIndex() && (samGzipIntFiles || mParams.outputParams().bam()) && noLongSequences && !delayMerge;
    final ReferenceRegions referenceRegions = mParams.outputParams().calibrateRegions();
    final OneShotTimer timer = new OneShotTimer("filterConcat");
    final int numIntermediateFiles = numberIntermediateFiles(tempFiles.length, mParams.numberThreads());
    final SimpleThreadPool pool = new SimpleThreadPool(Math.min(numIntermediateFiles, MAX_FILTERCONCAT_THREADS), mThreadNamePrefix + "FilterConcat", true);
    pool.enableBasicProgress(numIntermediateFiles);
    final File[] intermediate = new File[numIntermediateFiles];
    final File[] intermediateIndexes = createIndex ? new File[numIntermediateFiles] : null;
    final File[] intermediateCal = calibrate ? new File[numIntermediateFiles] : null;
    //    System.err.println("temp files: " + java.util.Arrays.toString(tempFiles));
    final Range[] regionRanges = groupRegions(tempFiles.length, numIntermediateFiles);

    for (int i = 0; i < regionRanges.length; ++i) {
      final File[] subFiles = new File[regionRanges[i].getLength()];
      System.arraycopy(tempFiles, regionRanges[i].getStart(), subFiles, 0, regionRanges[i].getLength());
      intermediate[i] = File.createTempFile("TEMP_FILTER_" + regionRanges[i].getStart() + "-" + regionRanges[i].getEnd() + "_", mParams.outputParams().bam() ? ".bam" : (".sam" + (samGzipIntFiles ? FileUtils.GZ_SUFFIX : "")), subFiles[0].getParentFile());
      //System.err.println(intermediate[i]);

      // Logic to prevent crash if we somehow get here with -Z option and BAM output
      // This should be prevented at the command line handling stage but wasn't in 3.3.2 or earlier
      final boolean intFilesCompress = mParams.outputParams().bam() || samGzipIntFiles;
      final OutputWrapper outWrapper = createStreams(numIntermediateFiles, intermediate, intermediateIndexes, intFilesCompress, createIndex, i);
      final OutputStream intCalStream;
      if (calibrate) {
        intermediateCal[i] = new File(intermediate[i].getParent(), intermediate[i].getName() + CommonFlags.RECALIBRATE_EXTENSION);
        intCalStream = FileUtils.createOutputStream(intermediateCal[i], false);
      } else {
        intCalStream = null;
      }
      final boolean writeHeader = i == 0;
      final AbstractSamResultsFilter filter = makeFilter();

      //TODO Should be put back when picard can read temp BAM file fragments (for svprep unmated input), remove equivalent lines from makeFilter implementations
      //See wiki RTG.Mapping for more information about workflow of the following.
      //filter.setBamOutput(mParams.outputParams().bam());
      filter.setHeader(header);
      filter.setWriteHeader(writeHeader);
      filter.setWriteBogusBamHeader(mParams.outputParams().unify() && mParams.outputParams().bam());

      //filter.setAddBamTerminator(mParams.outputParams().bam() && (i == numThreads - 1)); no longer required here, determined by outputstream
      pool.execute(new SubFilter(filter, outWrapper, intCalStream, referenceRegions, mParams.searchParams().reader().copy(), subFiles, header));
    }
    pool.terminate();

    postProcessIntermediateFiles(numIntermediateFiles, intermediate, intermediateIndexes, header, createIndex);

    if (!delayMerge) {
      final ArrayList<Long> dataFileSizes = new ArrayList<>();
      for (final File f : intermediate) {
        dataFileSizes.add(f.length());
      }
      final int tot = 1
              + ((samGzipIntFiles || mParams.outputParams().bam()) ? 1 : 0)
              + (calibrate ? 1 : 0);
      int cur = 0;
      final boolean keepTempFiles = GlobalFlags.isSet(CoreGlobalFlags.MAP_KEEP_TEMPORARY_FILES);
      Diagnostic.progress(mThreadNamePrefix + "ResultsConcat: Starting " + tot + " Jobs");
      FileUtils.catInSync(outFile, !keepTempFiles, intermediate);
      Diagnostic.progress(mThreadNamePrefix + "ResultsConcat: " + ++cur + "/" + tot + " Jobs Finished");
      if (createIndex) {
        mergeIndexes(mParams, outFile, intermediateIndexes, dataFileSizes);
        Diagnostic.progress(mThreadNamePrefix + "ResultsConcat: " + ++cur + "/" + tot + " Jobs Finished");
      }
      if (calibrate) {
        mergeCalibration(outFile, referenceRegions, intermediateCal);
        Diagnostic.progress(mThreadNamePrefix + "ResultsConcat: " + ++cur + "/" + tot + " Jobs Finished");
      }
    }
    timer.stopLog();
    return new FilterConcatIntermediateFiles(intermediate, intermediateCal, intermediateIndexes);
  }

  static int numberIntermediateFiles(int numberRegions, int numberThreads) {
    return Math.min(numberRegions, numberThreads);
  }

  static void mergeCalibration(File outFile, ReferenceRegions referenceRegions, File[] intermediateCal) throws IOException {
    final File outCal = new File(outFile.getPath() + CommonFlags.RECALIBRATE_EXTENSION);
    final Calibrator cal = new Calibrator(CovariateEnum.getCovariates(CovariateEnum.DEFAULT_COVARIATES, null), referenceRegions);
    for (final File f : intermediateCal) {
      cal.accumulate(f);
      if (!f.delete()) {
        Diagnostic.userLog("Failed to delete temporary file: " + f.getPath());
      }
    }
    cal.writeToFile(outCal);
  }

  static void mergeIndexes(NgsParams params, File outFile, File[] intermediateIndexes, List<Long> dataFileSizes) throws IOException {
    if (params.outputParams().bam()) {
      BamIndexMerge.mergeBamIndexFiles(indexFileName(outFile, true), Arrays.asList(intermediateIndexes), dataFileSizes);
    } else {
      TabixIndexMerge.mergeTabixFiles(indexFileName(outFile, false), Arrays.asList(intermediateIndexes), dataFileSizes);
    }
    for (final File f : intermediateIndexes) {
      if (f.exists() && !f.delete()) {
        Diagnostic.userLog("Failed to delete temporary file: " + f.getPath());
      }
    }
  }

  static Range[] groupRegions(int numberOfRegions, int numIntermediateFiles) {
    final Range[] regionRanges = new Range[numIntermediateFiles];
    for (int i = 0; i < numIntermediateFiles; ++i) {
      final int start = i * numberOfRegions / numIntermediateFiles;
      final int end = Math.min((i + 1) * numberOfRegions / numIntermediateFiles, numberOfRegions);
      regionRanges[i] = new Range(start, end);
    }
    return regionRanges;
  }

  /**
   * Return an <code>AbstractSamResultsFilter</code> suitable for filtering the temp files.
   *
   * @return an <code>AbstractSamResultsFilter</code> value
   */
  protected abstract AbstractSamResultsFilter makeFilter();


  @TestClass(value = {"com.rtg.ngs.TopNPairedEndOutputProcessorSyncTest"})
  static class SubFilter implements IORunnable {

    private final AbstractSamResultsFilter mFilter;
    private final OutputStream mOut;
    private final OutputStream mCalOut;
    private final SequencesReader mTemplate;
    private final File[] mFiles;
    private final SAMFileHeader mHeader;
    private final IndexingStreamCreator.IndexRunner mIndexRunner;
    private final ReferenceRegions mRegions;

    SubFilter(AbstractSamResultsFilter filter, OutputWrapper out, OutputStream calOut, ReferenceRegions referenceRegions, SequencesReader template, File[] files, SAMFileHeader header) {
      mFilter = filter;
      mOut = out.mOutputStream;
      mCalOut = calOut;
      mTemplate = template;
      mFiles = files;
      mHeader = header;
      mIndexRunner = out.mIndexRunner;
      mRegions = referenceRegions;
    }

    @Override
    public void run() throws IOException {
      final IORunnableProxy indexProxy = new IORunnableProxy(mIndexRunner);
      final Thread indexThread = new Thread(indexProxy);
      try (final OutputStream out = mOut;
           final OutputStream calOut = mCalOut;
           final SequencesReader template = mTemplate) {
        if (mIndexRunner != null) {
          indexThread.start();
        }
        mFilter.filterConcat(mHeader, out, calOut, mRegions, template, true, mFiles);
      } finally {
        if (mIndexRunner != null) {
          try {
            indexThread.join();
          } catch (final InterruptedException e) {
            throw new IOException("Execution was interrupted", e);
          } finally {
            indexProxy.checkError();
          }
        }
      }
      for (final File f : mFiles) {
        if (f.exists() && !f.delete()) {
          Diagnostic.userLog("Failed to delete temporary file: " + f.getPath());
        }
      }
    }
  }

  /**
   * Wrapper for an intermediate file output stream and an optional
   * index runner for the creation of indexes from the same output.
   */
  protected static class OutputWrapper {
    final OutputStream mOutputStream;
    final IndexingStreamCreator.IndexRunner mIndexRunner;

    /**
     * Constructor for wrapper.
     * @param outputStream the output stream for the intermediate file
     * @param indexRunner the index runner for producing the index for the intermediate file (can be null)
     */
    public OutputWrapper(OutputStream outputStream, IndexingStreamCreator.IndexRunner indexRunner) {
      mOutputStream = outputStream;
      mIndexRunner = indexRunner;
    }
  }

}
