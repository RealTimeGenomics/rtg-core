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

import static com.rtg.index.hash.ngs.ReadDecoder.PAIRED_END;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.ngs.blocking.ReadBlocker;
import com.rtg.ngs.blocking.ReadBlockerSync;
import com.rtg.ngs.tempstage.PairedTempFileWriterImpl;
import com.rtg.pairedend.ReadStatusListener;
import com.rtg.pairedend.SlidingWindowCollector;
import com.rtg.reader.NamesInterface;
import com.rtg.reader.PrereadType;
import com.rtg.sam.RecordIterator;
import com.rtg.sam.SamUtils;
import com.rtg.sam.SkipInvalidRecordsIterator;
import com.rtg.util.IORunnable;
import com.rtg.util.IORunnableProxy;
import com.rtg.util.Pair;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.diagnostic.Timer;
import com.rtg.util.io.FileUtils;
import com.rtg.variant.sv.ReadGroupStatsCalculator;
import com.rtg.variant.sv.UnmatedAugmenter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

/**
 * Thread safe version of super class
 */
public class TopNPairedEndOutputProcessorSync extends AbstractMapOutputProcessor {

  private static final boolean DUMB_N = false; //Boolean.valueOf(System.getProperty("rtg.dumbn", "false"));

  private final ReadBlocker mFreqBlockerLeft;
  private final ReadBlocker mFreqBlockerRight;

  //make dynamic
  private final List<Pair<HashingRegion, File>> mChildren;

  private int mChildCount = 0;

  private final boolean mOutputUnmated;
  private final boolean mOutputUnmapped;

  private UptoNStore mTopN;

  /**
   * Construct an output processor
   * @param param parameters
   * @param stats map to put statistics into
   * @param outputUnmated true if unmated/unmapped should be output
   * @param outputUnmapped  true if unmapped should be output
   * @throws IOException if sam has problems setting up.
   */
  public TopNPairedEndOutputProcessorSync(NgsParams param, MapStatistics stats, boolean outputUnmated, boolean outputUnmapped) throws IOException {
    super(param, stats, true, outputUnmated);
    final long sequences = param.buildFirstParams().numberSequences();

    // These blockers are for counting hits per read per side, for high frequency filtering purposes
    mFreqBlockerLeft = new ReadBlockerSync(sequences, param.readFreqThreshold(), "left hits");
    mFreqBlockerRight = new ReadBlockerSync(sequences, param.readFreqThreshold(), "right hits");
    mChildren = new ArrayList<>();
    mOutputUnmated = outputUnmated;
    mOutputUnmapped = outputUnmapped;
    if (mOutputUnmated) {
      final UptoNStore nStore;
      if (DUMB_N) {
        Diagnostic.developerLog("Using DumbN: " + mParams.outputParams().filter().topN());
        nStore = new DeduplicatingNStore(sequences * 2, param.searchParams().reader().numberSequences(), mParams.outputParams().filter().topN(), param.searchParams().reader().maxLength(), Math.max(param.buildFirstParams().maxLength(), param.buildSecondParams().maxLength()));
      } else {
        Diagnostic.developerLog("Using TopN: " + mParams.outputParams().filter().topN());
        nStore = new TopNImplementation(sequences * 2, param.searchParams().reader().numberSequences(), mParams.outputParams().filter().topN(), param.searchParams().reader().maxLength(), Math.max(param.buildFirstParams().maxLength(), param.buildSecondParams().maxLength()));
      }
      mTopN = new UptoNStoreSync(nStore);
    } else {
      mTopN = null;
    }
    //mChildren = new PairedEndOutputProcessor[mParams.numberThreads()];
  }

  private static synchronized void createDir(File dir) throws IOException {
    if (!dir.exists() && !dir.mkdirs()) {
      throw new IOException("unable to create directory: " + dir);
    }
  }

  @Override
  public synchronized OutputProcessor threadClone(HashingRegion region) throws IOException {
    super.threadClone(region);
    final int currentChild = mChildCount++;

    // Check that output directory exists
    final File dir = mParams.outputParams().directory();
    createDir(dir);
    final File out = TopNPairedEndOutputProcessorSync.determineTempFile(mParams, currentChild);
    final OutputStream outStream = FileUtils.createOutputStream(out, true, false);

    final PairedTempFileWriterImpl sam = new PairedTempFileWriterImpl(mParams,  mUnmappedTracker, mSharedResources);
    sam.initialiseMated(outStream);
    sam.setClipRegion(region);
    final SlidingWindowCollector collector = new SlidingWindowCollector(mParams.maxFragmentLength(), mParams.minFragmentLength(), mParams.pairOrientation(), sam,
        mSharedResources, mParams.outputParams().calibrateRegions());
    final InnerPairedEndOutputProcessor ipeop = new InnerPairedEndOutputProcessor(sam, collector, mFreqBlockerLeft, mFreqBlockerRight, mOutputUnmated ? mTopN : null, mUnmappedTracker, region);
    mChildren.add(new Pair<>(ipeop.getRegion(), out));
    Diagnostic.developerLog("InnerPairedEndOutputProcessor:" + ipeop + " region:" + region);
    return ipeop;
  }

  @Override
  public void threadFinish() {
    throw new UnsupportedOperationException("not allowed for this class");
  }

  static class RegionFileComparator implements Comparator<Pair<HashingRegion, File>>, Serializable {
    @Override
    public int compare(Pair<HashingRegion, File> first, Pair<HashingRegion, File> second) {
      return first.getA().compareTo(second.getA());
    }
  }

  @Override
  public void finish() throws IOException {
    sortRegions();
    if (mTopN != null) {
      Diagnostic.developerLog(mTopN.histogram());
    }
    Collections.sort(mChildren, new RegionFileComparator());
    final File[] outputFiles = new File[mChildren.size()];
    for (int i = 0; i < outputFiles.length; ++i) {
      outputFiles[i] = mChildren.get(i).getB();
    }
    // Filter and concatenate the intermediate files. All thresholds etc are now known.
    final File outFile;
    if (mParams.outputParams().bam()) {
      outFile = mParams.outputParams().resultStreamHandler().file(NgsOutputParams.MATED_BAM_FILE_NAME);
    } else {
      outFile = mParams.outputParams().isCompressOutput()
      ? mParams.outputParams().resultStreamHandler().file(NgsOutputParams.MATED_SAM_FILE_NAME + FileUtils.GZ_SUFFIX)
      : mParams.outputParams().resultStreamHandler().file(NgsOutputParams.MATED_SAM_FILE_NAME);
    }
    preProcessMatedStatistics(mUnmappedTracker, mSharedResources.getBlocker()) ;

    final FilterConcatIntermediateFiles unmatedIntFiles;
    if (mOutputUnmated) {
      unmatedIntFiles = writeUnmated();
    } else {
      unmatedIntFiles = null;
    }

    final FilterConcatIntermediateFiles matedIntFiles = new MatedMulticoreFilterConcat(mParams, mSharedResources.getBlocker(), mFreqBlockerLeft, mFreqBlockerRight, mSharedResources.names(), mParams.useTopRandom() ? mSharedResources.getPairedEndTopRandom().getRecords() : null, mStatsMerger, mReportMerger)
      .filterConcat(outputFiles, outFile, mSharedResources.getHeader(), mParams.outputParams());
    if (mSharedResources.getPairedEndTopRandom() != null) {
      mSharedResources.getPairedEndTopRandom().finish();
    }

    final FilterConcatIntermediateFiles unmappedIntFiles;
    if (mOutputUnmapped) {
      unmappedIntFiles = writeUnmapped(!mParams.outputParams().unify(), false, false);
    } else {
      unmappedIntFiles = null;
    }
    if (mParams.outputParams().unify()) {
      whizBangUnify(unmappedIntFiles, matedIntFiles, unmatedIntFiles);
      //unifyAlignmentOutput(alignmentFiles);
    }

    if (mStatsMerger != null) {
      try (OutputStream rgOut = FileUtils.createOutputStream(mParams.outputParams().resultStreamHandler().file(UnmatedAugmenter.DEFAULT_RGSTATS_FILENAME), false)) {
        mStatsMerger.blend().writeReadGroupStats(rgOut);
      }
    }
    mUnmappedTracker.calculateStatistics(true, false);
    mReportMerger.blendReportData().write(new File(mParams.outputParams().directory(), MapReportData.MAP_REPORT_FILE_NAME));
  }

  protected FilterConcatIntermediateFiles writeUnmated() throws IOException {
    Diagnostic.userLog("Extracting hits for unmated reads");
    Diagnostic.progress("UnmatedInit: Starting 1 Jobs");
    final Timer unmatedOutputTimer = new Timer("UnmatedOutput");
    unmatedOutputTimer.start();
    final MatchResult results = TopNPairedEndOutputProcessorSync.organizeUnmatedResults(mUnmappedTracker.mReadIdStatus, mTopN);

    mTopN = null;
    final boolean paired = true;
    final HashingRegion[] regions = new HashingRegion[mChildren.size()];
    for (int i = 0; i < mChildren.size(); ++i) {
      regions[i] = mChildren.get(i).getA();
    }
    final FilterConcatIntermediateFiles unmatedIntFiles = getNonMatedFilterConcatIntermediateFiles(results, paired, regions);
    unmatedOutputTimer.stop();
    unmatedOutputTimer.log();
    return unmatedIntFiles;
  }

  @Override
  protected FilterConcatIntermediateFiles filterConcatNonMated(MapQScoringReadBlocker blockerLeft, MapQScoringReadBlocker blockerRight, File[] tempFiles, SingleEndTopRandomImplementation.HitRecord[] hitsToKeep, NamesInterface templateNames, File outFile) throws IOException {
    preProcessUnMatedStatistics(mUnmappedTracker, blockerLeft, blockerRight, mFreqBlockerLeft, mFreqBlockerRight);
    return new UnmatedMulticoreFilterConcat(mParams, blockerLeft, blockerRight, mFreqBlockerLeft, mFreqBlockerRight, hitsToKeep, templateNames, mAugmenterMerger, mStatsMerger, mReportMerger)
            .filterConcat(tempFiles, outFile, mSharedResources.getHeader(), mParams.outputParams());
  }

  private void preProcessUnMatedStatistics(ReadStatusTracker unmappedTracker, MapQScoringReadBlocker asBlockerLeft, MapQScoringReadBlocker asBlockerRight, ReadBlocker freqBlockerLeft, ReadBlocker freqBlockerRight) {
    Diagnostic.progress("UnmatedPreprocess: Starting 1 Jobs");
    for (int i = 0; i < mParams.buildFirstParams().numberSequences(); ++i) {
      if (unmappedTracker.getStatus(i, ReadStatusTracker.MATED_FIRST) && unmappedTracker.getStatus(i, ReadStatusTracker.MATED_SECOND)) {
        continue;
      }
      if (!asBlockerLeft.isBlocked(i) && !freqBlockerLeft.isBlocked(i)) {
        final int count = asBlockerLeft.getCount1(i);
        if (count > 0) {
          if (count == 1) {
            unmappedTracker.addStatus(i, ReadStatusTracker.UNIQUELY_MAPPED_FIRST);
          }
          unmappedTracker.addStatus(i, ReadStatusTracker.UNMATED_FIRST);
        }
      }
      if (!asBlockerRight.isBlocked(i) && !freqBlockerRight.isBlocked(i)) {
        final int count = asBlockerRight.getCount1(i);
        if (count > 0) {
          if (count == 1) {
            unmappedTracker.addStatus(i, ReadStatusTracker.UNIQUELY_MAPPED_SECOND);
          }
          unmappedTracker.addStatus(i, ReadStatusTracker.UNMATED_SECOND);
        }
      }
    }
    Diagnostic.progress("UnmatedPreprocess: 1/1 Jobs Finished");

  }

  private void preProcessMatedStatistics(ReadStatusTracker unmappedTracker, MapQScoringReadBlocker blocker) {
    Diagnostic.progress("MatedPreprocess: Starting 1 Jobs");
    for (int i = 0; i < mParams.buildFirstParams().numberSequences(); ++i) {
      if (blocker.isBlocked(i) || mFreqBlockerLeft.isBlocked(i) || mFreqBlockerRight.isBlocked(i)) {
        //System.err.println("skipping ... " + i);
        continue;
      }
      final int count = blocker.getCount1(i);
      if (count > 0) {
        if (count == 1) {
          unmappedTracker.addStatus(i, ReadStatusTracker.UNIQUELY_MAPPED_FIRST);
          unmappedTracker.addStatus(i, ReadStatusTracker.UNIQUELY_MAPPED_SECOND);
        }
        unmappedTracker.addStatus(i, ReadStatusTracker.MATED_FIRST);
        unmappedTracker.addStatus(i, ReadStatusTracker.MATED_SECOND);
      }
    }
    Diagnostic.progress("MatedPreprocess: 1/1 Jobs Finished");
  }

  static class MatedMulticoreFilterConcat extends AbstractMulticoreFilterConcat {
    final MapQScoringReadBlocker mXaBlocker;
    final ReadBlocker mFreqBlockerLeft;
    final ReadBlocker mFreqBlockerRight;
    final NamesInterface mTemplateNames;
    final PairedTopRandomImplementation.HitRecord[] mHitsToKeep;
    final ReadGroupStatsCalculator.Merger mStatsMerger;
    final MapReportData.Merger mReportMerger;

    MatedMulticoreFilterConcat(NgsParams params, MapQScoringReadBlocker xaBlocker, ReadBlocker freqBlockerLeft, ReadBlocker freqBlockerRight, NamesInterface templateNames, PairedTopRandomImplementation.HitRecord[] hitsToKeep, ReadGroupStatsCalculator.Merger statsMerger, MapReportData.Merger reportMerger) {
      super(params, "Mated");
      mXaBlocker = new MapQScoringReadBlocker(xaBlocker);    // Read-only from here, so no need for synchronization
      mFreqBlockerLeft = new ReadBlocker(freqBlockerLeft);   // Read-only from here, so no need for synchronization
      mFreqBlockerRight = new ReadBlocker(freqBlockerRight); // Read-only from here, so no need for synchronization
      mTemplateNames = templateNames;
      mHitsToKeep = hitsToKeep;
      mStatsMerger = statsMerger;
      mReportMerger = reportMerger;
    }

    @Override
    protected AbstractSamResultsFilter makeFilter() {
      final MatedSamResultsFilter filter;
      final long readIdOffset = Math.max(0, mParams.buildFirstParams().readerRestriction().getStart());
      final String readGroupId = mParams.outputParams().readGroup() != null ? mParams.outputParams().readGroup().getReadGroupId() : null;
      filter = new MatedSamResultsFilter(mXaBlocker, mFreqBlockerLeft, mFreqBlockerRight, mParams.buildFirstParams().reader().copy(), mParams.buildSecondParams().reader().copy(), mParams.buildFirstParams().reader().getPrereadType() == PrereadType.CG, readIdOffset, readGroupId, mParams.legacyCigars());
      filter.setBamOutput(mParams.outputParams().bam());
      filter.setStatsCalculator(mStatsMerger);
      filter.setMapReportData(mReportMerger);
      if (mParams.useTopRandom()) {
        filter.setTemplateNames(mTemplateNames);
        filter.setHitsToKeep(mHitsToKeep);
      }
      if (mParams.outputParams().outputReadNames()) {
        try {
          filter.setReadNames(mParams.buildFirstParams().reader().names());
        } catch (final IOException e) {
          filter.setReadNames(null);
          Diagnostic.warning("Failed to retrieve read names, using read id instead.");
        }
      }
      return filter;
    }
  }

  static class UnmatedMulticoreFilterConcat extends AbstractMulticoreFilterConcat {
    final MapQScoringReadBlocker mAsBlockerLeft;
    final MapQScoringReadBlocker mAsBlockerRight;
    final ReadBlocker mFreqBlockerLeft;
    final ReadBlocker mFreqBlockerRight;
    final SingleEndTopRandomImplementation.HitRecord[] mHitsToKeep;
    final NamesInterface mTemplateNames;
    final UnmatedAugmenter.Merger mAugmenterMerger;
    final ReadGroupStatsCalculator.Merger mStatsMerger;
    final MapReportData.Merger mReportMerger;

    UnmatedMulticoreFilterConcat(NgsParams params, MapQScoringReadBlocker asBlockerLeft, MapQScoringReadBlocker asBlockerRight, ReadBlocker freqBlockerLeft, ReadBlocker freqBlockerRight, SingleEndTopRandomImplementation.HitRecord[] hitsToKeep, NamesInterface templateNames, UnmatedAugmenter.Merger augmenterMerger, ReadGroupStatsCalculator.Merger statsMerger, MapReportData.Merger reportMerger) {
      super(params, "Unmated");
      mAsBlockerLeft = new MapQScoringReadBlocker(asBlockerLeft);    // Read-only from here, so no need for synchronization
      mAsBlockerRight = new MapQScoringReadBlocker(asBlockerRight);  // Read-only from here, so no need for synchronization
      mFreqBlockerLeft = new ReadBlocker(freqBlockerLeft);           // Read-only from here, so no need for synchronization
      mFreqBlockerRight = new ReadBlocker(freqBlockerRight);         // Read-only from here, so no need for synchronization
      mHitsToKeep = hitsToKeep;
      mTemplateNames = templateNames;
      mAugmenterMerger = augmenterMerger;
      mStatsMerger = statsMerger;
      mReportMerger = reportMerger;
    }
    @Override
    protected AbstractSamResultsFilter makeFilter() {
      final long readIdOffset = Math.max(0, mParams.buildFirstParams().readerRestriction().getStart());
      final String readGroupId = mParams.outputParams().readGroup() != null ? mParams.outputParams().readGroup().getReadGroupId() : null;
      final UnmatedSamResultsFilter filter = new UnmatedSamResultsFilter(mAsBlockerLeft, mAsBlockerRight, mFreqBlockerLeft, mFreqBlockerRight, readIdOffset, mParams.buildFirstParams().reader().copy(), mParams.buildSecondParams().reader().copy(), readGroupId, mParams.legacyCigars(), mAugmenterMerger);
      filter.setBamOutput(mAugmenterMerger == null && mParams.outputParams().bam()); //TODO replace hack to deal with the inability of picard to read in the temp BAM file fragments properly by using SAM on initial output with actually using BAM when picard can handle it
      filter.setReadersAndType(mParams.buildFirstParams().reader().copy(), mParams.buildSecondParams().reader().copy(), mParams.buildFirstParams().reader().getPrereadType() == PrereadType.CG);
      filter.setStatsCalculator(mStatsMerger);
      filter.setMapReportData(mReportMerger);
      if (mParams.useTopRandom()) {
        filter.setHitsToKeep(mHitsToKeep);
        filter.setTemplateNames(mTemplateNames);
      }
      if (mParams.outputParams().outputReadNames()) {
        try {
          filter.setReadNames(mParams.buildFirstParams().reader().names());
        } catch (final IOException e) {
          filter.setReadNames(null);
          Diagnostic.warning("Failed to retrieve read names, using read id instead.");
        }
      }
      return filter;
    }

    @Override
    protected OutputWrapper createStreams(int numThreads, File[] intermediate, File[] intermediateIndexes, boolean samGzipIntFiles, boolean createIndex, int i) throws IOException {
      //Don't want to index files yet due to need to post-process intermediate files.
      if (mAugmenterMerger != null) {
        return new OutputWrapper(FileUtils.createOutputStream(intermediate[i], samGzipIntFiles, false, samGzipIntFiles && (i == numThreads - 1)), null);
      }
      return super.createStreams(numThreads, intermediate, intermediateIndexes, samGzipIntFiles, createIndex, i);
    }

    @Override
    protected void postProcessIntermediateFiles(int numThreads, File[] intermediate, File[] intermediateIndexes, SAMFileHeader header, boolean createIndex) throws IOException {
      if (mAugmenterMerger == null) {
        return;
      }
      //read in and update with pseudo-pairing information, and re-write intermediate files, write index files if needed also
      final SimpleThreadPool pool = new SimpleThreadPool(Math.min(numThreads, MAX_FILTERCONCAT_THREADS), mThreadNamePrefix + "FilterConcat-PostProcessing", true);
      pool.enableBasicProgress(numThreads);
      final File[] intermediateTempFiles = new File[numThreads];
      final boolean samGzipIntFiles = mParams.outputParams().isCompressOutput();
      final UnmatedAugmenter au = mAugmenterMerger.blend();
      for (int i = 0; i < numThreads; ++i) {
        intermediateTempFiles[i] = new File(intermediate[i].getParent(), intermediate[i].getName() + ".augment");
        final OutputWrapper outWrapper = super.createStreams(numThreads, intermediateTempFiles, intermediateIndexes, samGzipIntFiles, createIndex, i);
        final boolean writeHeader = i == 0;
        final boolean bogoHeader = mParams.outputParams().bam() && mParams.outputParams().unify();
        pool.execute(new SubAugment(outWrapper, intermediate[i], intermediateTempFiles[i], au, header, writeHeader, mParams.outputParams().bam(), bogoHeader));
      }
      pool.terminate();
    }

    static class SubAugment implements IORunnable {
      private final OutputStream mOut;
      private final IORunnableProxy mIndexProxy;
      private final File mIntermediateFile;
      private final File mIntermediateAugmentFile;
      private final UnmatedAugmenter mAugmenter;
      private final SAMFileHeader mHeader;
      private final boolean mWriteHeader;
      private final boolean mBamOutput;
      private final boolean mWriteBogusHeader;

      /**
       * @param out output wrapper containing output stream to write to
       * @param intermediate the intermediate file to read from. May have a missing or fake header.
       * @param intermediateAugment the output file we're writing to in <code>out</code>, so we can move this file over <code>intermediate</code> after augmenting
       * @param augmenter the augmenter instance to use
       * @param header the real sam header. Used to override the header while reading, which may be missing or fake.
       * @param writeHeader true to write a header to this file. Omitting the header is desirable for intermediate files after the first. Ignored if <code>writeBogusHeader</code> is true.
       * @param bam true to use bam writers for output
       * @param writeBogusHeader true to write out minimal headers instead of the normal header. Usually for in unified BAM output case (Picard will later read these and need a header to detect BAM file type)
       */
      SubAugment(OutputWrapper out, File intermediate, File intermediateAugment, UnmatedAugmenter augmenter, SAMFileHeader header, boolean writeHeader, boolean bam, boolean writeBogusHeader) {
        mOut = out.mOutputStream;
        mIndexProxy = out.mIndexRunner == null ? null : new IORunnableProxy(out.mIndexRunner);
        mIntermediateFile = intermediate;
        mIntermediateAugmentFile = intermediateAugment;
        mAugmenter = augmenter;
        mHeader = header;
        mWriteHeader = writeHeader;
        mBamOutput = bam;
        mWriteBogusHeader = writeBogusHeader;
      }

      @Override
      public void run() throws IOException {
        final Thread indexThread = mIndexProxy == null ? null : new Thread(mIndexProxy);
        try (final OutputStream out = mOut) {
          if (indexThread != null) {
            indexThread.start();
          }

          final SAMFileHeader basicHeader = new SAMFileHeader();
          basicHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
          try (SamReader reader = SamUtils.makeSamReader(FileUtils.createFileInputStream(mIntermediateFile, false), null, mHeader)) {
            try (SAMFileWriter writer = getSAMFileWriter(basicHeader, out)) {
              final RecordIterator<SAMRecord> it = new SkipInvalidRecordsIterator(mIntermediateFile.getPath(), reader, true);
              while (it.hasNext()) {
                final SAMRecord rec = it.next();
                mAugmenter.updateUnmatedRecord(rec);
                writer.addAlignment(rec);
              }
            }
          }
        } finally {
          if (indexThread != null) {
            try {
              indexThread.join();
            } catch (final InterruptedException e) {
              throw new IOException("Execution was interrupted", e);
            } finally {
              mIndexProxy.checkError();
            }
          }
        }
        if (!(mIntermediateFile.canWrite() && mIntermediateFile.delete() && mIntermediateAugmentFile.renameTo(mIntermediateFile))) {
          throw new SlimException("Unable to rename intermediate file: \"" + mIntermediateAugmentFile + "\" to \"" + mIntermediateFile + "\"");
        }
      }

      private SAMFileWriter getSAMFileWriter(SAMFileHeader basicHeader, OutputStream out) {
        final SAMFileWriter writer;
        if (mBamOutput && mWriteBogusHeader) {
          writer = new SAMFileWriterFactory().makeBAMWriter(basicHeader, true, out, true, false, true);
        } else if (mBamOutput) {
          writer = new SAMFileWriterFactory().makeBAMWriter(mHeader, true, out, mWriteHeader, false, true);
        } else {
          writer = new SAMFileWriterFactory().makeSAMWriter(mHeader, true, out, mWriteHeader);
        }
        return writer;
      }
    }
  }

  @Override
  public void close() throws IOException {
    mFreqBlockerLeft.close();
    mFreqBlockerRight.close();
    mSharedResources.close();
  }


  @Override
  public void process(long templateId, String frame, int readId, int tStart, int score, int scoreIndel) {
    throw new UnsupportedOperationException("Calls should be made on thread local clones");
  }

  /**
   * Creates a temporary file for use with unfiltered mated sam output.
   * @param param contains settings used to locate the temporary directory and filename extension
   * @param currentChild a numeric indicator used as part of the filename
   * @return the temp file
   * @throws IOException if an IO error occurs.
   */
  public static synchronized File determineTempFile(final NgsParams param, final int currentChild) throws IOException {
    final File tempDir = param.outputParams().tempFilesDirectory();
    // check that the temporary output directory exists
    if (!tempDir.exists() && !tempDir.mkdirs()) {
      throw new IOException("Could not create temporary directory: " + tempDir.getPath());
    }
    final boolean gzipTempFiles = param.outputParams().isCompressOutput();
    return param.outputParams().resultStreamHandler().tempFile(AbstractMapOutputProcessor.TEMP_SAM_ALIGNMENT_NAME + currentChild
            + (gzipTempFiles ? FileUtils.GZ_SUFFIX : ""));
  }

  private static MatchResult organizeUnmatedResults(final int[] readIdStatus, final UptoNStore uptoN) {
    int guessCount = 0;
    for (final int readStatus : readIdStatus) {
      //System.err.println(i + " READ ID STATUS = " + mReadIdStatus[i]);
      if (((readStatus & ReadStatusTracker.MATED_FIRST) == 0) && ((readStatus & ReadStatusTracker.MATED_SECOND) == 0)) {
        guessCount += 2;
      }
    }
    final MatchResult results = new MatchResult(guessCount);

    for (int i = 0; i < readIdStatus.length; ++i) {
      //System.err.println(i + " READ ID STATUS = " + mReadIdStatus[i]);
      if (((readIdStatus[i] & ReadStatusTracker.MATED_FIRST) == 0) && ((readIdStatus[i] & ReadStatusTracker.MATED_SECOND) == 0)) {
        uptoN.setResults(results, i * 2);
        uptoN.setResults(results, i * 2 + 1);
      }
    }
    results.sort();
    return results;
  }

  private static class InnerPairedEndOutputProcessor implements OutputProcessor {
    private final boolean mOutputUnmated;
    private final UptoNStore mUptoN;
    private final ReadStatusListener mListener;
    private final HashingRegion mRegion;
    private final PairedEndOutputProcessor mPairedEndOutputProcessor;

    protected final ReadBlocker mFreqBlockerLeft;
    protected final ReadBlocker mFreqBlockerRight;
    private int mProcessCalls = 0;

    InnerPairedEndOutputProcessor(PairedTempFileWriterImpl sam, SlidingWindowCollector collector, ReadBlocker left, ReadBlocker right, UptoNStore topn, ReadStatusListener listener, HashingRegion region) {
      mPairedEndOutputProcessor = new PairedEndOutputProcessor(sam, collector);
      mUptoN = topn;
      mOutputUnmated = mUptoN != null;
      mListener = listener;
      mRegion = region;
      mFreqBlockerLeft = left;
      mFreqBlockerRight = right;
    }

    /**
     * @return the HashingRegion this output processor is processing
     */
    public HashingRegion getRegion() {
      return mRegion;
    }

    /** Not supported in thread local version */
    @Override
    public OutputProcessor threadClone(final HashingRegion region) {
      throw new UnsupportedOperationException("Not supported.");
    }

    @Override
    public void threadFinish() throws IOException {
      try {
        finish();
      } finally {
        close();
      }
    }

    @Override
    public void close() throws IOException {
      mPairedEndOutputProcessor.close();
    }

    @Override
    public void finish() throws IOException {
      mPairedEndOutputProcessor.finish();
      Diagnostic.developerLog("Child finish TopNPEOPChild process calls: " + mProcessCalls);
    }

    @Override
    public void process(long templateId, String frame, int readId, int tStart, int score, int scoreIndel) throws IOException {
      final int dReadId = PAIRED_END.decode(readId);
      final boolean isFirst = PAIRED_END.isFirst(readId);
      if (isFirst) {
        if (mFreqBlockerLeft.isBlocked(dReadId)) {
          mListener.addStatus(dReadId, ReadStatusTracker.BLOCKED_FIRST);
          return;
        }
      } else {
        if (mFreqBlockerRight.isBlocked(dReadId)) {
          mListener.addStatus(dReadId, ReadStatusTracker.BLOCKED_SECOND);
          return;
        }
      }
      ++mProcessCalls;
      final boolean reverse = frame.startsWith("R");
      mListener.addStatus(dReadId, isFirst ? ReadStatusTracker.MATCHED_FIRST : ReadStatusTracker.MATCHED_SECOND);
      mPairedEndOutputProcessor.process(templateId, reverse, readId, tStart);
      if (mRegion.isInRange(templateId, tStart)) {
        if (isFirst) {
          mFreqBlockerLeft.increment(dReadId);
          if (mFreqBlockerLeft.isBlocked(dReadId)) {
            mListener.addStatus(dReadId, ReadStatusTracker.BLOCKED_FIRST);
            return;
          }
        } else {
          mFreqBlockerRight.increment(dReadId);
          if (mFreqBlockerRight.isBlocked(dReadId)) {
            mListener.addStatus(dReadId, ReadStatusTracker.BLOCKED_SECOND);
            return;
          }
        }
        if (mOutputUnmated) {
            mUptoN.process(templateId, reverse, readId, tStart, scoreIndel);
        }
      }
    }
  }

}
