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
package com.rtg.sam;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collection;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.Recalibrate;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.tabix.IndexingStreamCreator;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.BlockCompressedOutputStream;

/**
 */
@TestClass("com.rtg.sam.SamMergeCliTest")
public class SamMerger {

  private final boolean mCreateIndex;
  private final boolean mGzip;
  private final boolean mLegacy;
  private final int mNumberThreads;
  private final SamFilterParams mFilterParams;
  private final boolean mAddProgramRecord;
  private final boolean mDeleteInputFiles;

  /**
   * @param createIndex true to create index for output files
   * @param gzip true to zip output files
   * @param legacy output cigars in legacy format
   * @param numberThreads number of threads to use for merging
   * @param filterParams filter options
   * @param addProgramRecord true if this merge is being run from command line and should add a program record to the output SAM/BAM header
   * @param deleteInputFiles true if should delete input SAM/BAM files and associated index and calibration files when merging
   */
  public SamMerger(boolean createIndex, boolean gzip, boolean legacy, int numberThreads, SamFilterParams filterParams, boolean addProgramRecord, boolean deleteInputFiles) {
    this.mCreateIndex = createIndex;
    this.mGzip = gzip;
    this.mLegacy = legacy;
    this.mNumberThreads = numberThreads;
    this.mFilterParams = filterParams;
    this.mAddProgramRecord = addProgramRecord;
    this.mDeleteInputFiles = deleteInputFiles;
  }

  /**
   * Merges given SAM and/or BAM files into one using parameters given to this merger.
   *
   * @param samFiles files to merge
   * @param calibrationFiles corresponding calibration files to merge (may be empty)
   * @param output output file, null to use supplied output stream
   * @param out output stream to write sam to (if {@code output} is null), otherwise statistics are written here. may be null
   * @param header use this header instead of any found in the SAM/BAM file (useful if file has no header, or needs changing).
   * @param writeHeader true if SAM/BAM header should be written to output file
   * @param terminateBlockedGzip true if BAM or block compressed SAM should have terminator block
   * @throws java.io.IOException if an IO error occurs
   */
  public void mergeSamFiles(Collection<File> samFiles, Collection<File> calibrationFiles, File output, OutputStream out, SAMFileHeader header, boolean writeHeader, boolean terminateBlockedGzip) throws IOException {
    if (output != null && calibrationFiles.size() > 0 && calibrationFiles.size() != samFiles.size()) {
      Diagnostic.warning("Number of calibration files does not match number of SAM files, will not merge calibration files.");
    }
    final boolean outputBam = (output != null) && output.getName().endsWith(SamUtils.BAM_SUFFIX);
    final File alignmentOutputFile = output != null ? (outputBam ? output : SamUtils.getZippedSamFileName(mGzip, output)) : null;
    final TabixIndexer.IndexerFactory indexerFactory = outputBam ? null : new TabixIndexer.SamIndexerFactory();
    final long recordsIn;
    final long recordsOut;
    final SingletonPopulatorFactory<SAMRecord> pf = new SingletonPopulatorFactory<>(new SamRecordPopulator());
    try (final ThreadedMultifileIterator<SAMRecord> it = new ThreadedMultifileIterator<>(samFiles, mNumberThreads, pf, mFilterParams, header)) {
      header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
      if (mAddProgramRecord) {
        SamUtils.addProgramRecord(header);
      }
      SamUtils.updateRunId(header);
      try (IndexingStreamCreator streamHandler = new IndexingStreamCreator(alignmentOutputFile, out, mGzip, indexerFactory, mCreateIndex)) {
        try (SAMFileWriter writer = getSAMFileWriter(streamHandler, header, outputBam, writeHeader, terminateBlockedGzip)) {
          while (it.hasNext()) {
            final SAMRecord rec = it.next();
            if (mLegacy) {
              SamUtils.convertToLegacyCigar(rec);
            }
            writer.addAlignment(rec);
          }
        }
      }
      recordsIn = it.getTotalRecordsCount();
      recordsOut = it.getTotalRecordsCount() - it.getFilteredRecordsCount() - it.getDuplicateRecordsCount() - it.getInvalidRecordsCount();
    }
    if (output != null) {
      if (calibrationFiles.size() > 0 && calibrationFiles.size() == samFiles.size()) {
        final Calibrator c = new Calibrator(Calibrator.getCovariateSet(calibrationFiles.iterator().next()), null);
        for (final File f : calibrationFiles) {
          c.accumulate(f);
        }
        final File actualOut = outputBam ? output : SamUtils.getZippedSamFileName(mGzip, output);
        c.writeToFile(new File(actualOut.getParent(), actualOut.getName() + Recalibrate.EXTENSION));
      }
      if (out != null) {
        out.write(("SAM records read:    " + recordsIn + StringUtils.LS).getBytes());
        out.write(("SAM records written: " + recordsOut + StringUtils.LS).getBytes());
      }
      if (mDeleteInputFiles) {
        for (final File f : samFiles) {
          final File indexName;
          if (f.getName().endsWith(SamUtils.BAM_SUFFIX)) {
            indexName = BamIndexer.indexFileName(f);
          } else {
            indexName = TabixIndexer.indexFileName(f);
          }
          if (indexName.exists() && !indexName.delete()) {
            Diagnostic.warning("Failed to delete file: " + indexName.getPath());
          }
          if (!f.delete()) {
            Diagnostic.warning("Failed to delete file: " + f.getPath());
          }
        }
        for (final File f : calibrationFiles) {
          if (f.exists() && !f.delete()) {
            Diagnostic.warning("Failed to delete file: " + f.getPath());
          }
        }
      }
    }
  }

  private static final int GZIP_LEVEL = GlobalFlags.getIntegerValue(CoreGlobalFlags.GZIP_LEVEL);

  /**
   * Get a SAM/BAM writer potentially with indexing on the fly
   * @param streamHandler stream handler
   * @param header header to output
   * @param outputBam if true BAM otherwise SAM
   * @param writeHeader should the header be written
   * @param terminateBlockGzip should the terminator be written
   * @return writer
   * @throws IOException if an I/O error occurs.
   */
  public static SAMFileWriter getSAMFileWriter(IndexingStreamCreator streamHandler, SAMFileHeader header, boolean outputBam, boolean writeHeader, boolean terminateBlockGzip) throws IOException {
    final OutputStream intStream = streamHandler.createStreamsAndStartThreads(header != null ? header.getSequenceDictionary().size() : -1, writeHeader, terminateBlockGzip);
    final SAMFileWriter writer;
    if (outputBam) {
      writer = new SAMFileWriterFactory().makeBAMWriter(header, true, new BlockCompressedOutputStream(intStream, null, GZIP_LEVEL, terminateBlockGzip), writeHeader, false /* ignored */, true);
    } else {
      writer = new SAMFileWriterFactory().makeSAMWriter(header, true, intStream, writeHeader);
    }
    return writer;
  }

}
