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
import com.rtg.launcher.CommonFlags;
import com.rtg.reader.SequencesReader;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

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
  private boolean mRenameWithRg;

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

  void setRenameWithRg(boolean rename) {
    mRenameWithRg = rename;
  }

  /**
   * Merges given SAM and/or BAM files into one using parameters given to this merger.
   *
   * @param samFiles files to merge
   * @param calibrationFiles corresponding calibration files to merge (may be empty)
   * @param output output file, or "-" to use supplied output stream
   * @param out output stream to write sam to (if {@code output} is "-"), otherwise statistics are written here. may be null
   * @param reference must be non-null for CRAM support
   * @param header use this header instead of any found in the SAM/BAM file (useful if file has no header, or needs changing).
   * @param writeHeader true if SAM/BAM header should be written to output file
   * @param terminateBlockedGzip true if BAM or block compressed SAM should have terminator block
   * @throws java.io.IOException if an IO error occurs
   */
  public void mergeSamFiles(Collection<File> samFiles, Collection<File> calibrationFiles, File output, OutputStream out, SequencesReader reference, SAMFileHeader header, boolean writeHeader, boolean terminateBlockedGzip) throws IOException {
    mergeSamFiles(samFiles, calibrationFiles, output, out, reference, header, writeHeader ? header : null, terminateBlockedGzip);
  }

  /**
   * Merges given SAM and/or BAM files into one using parameters given to this merger.
   *
   * @param samFiles files to merge
   * @param calibrationFiles corresponding calibration files to merge (may be empty)
   * @param output output file, or "-" to use supplied output stream
   * @param out output stream to write sam to (if {@code output} is "-"), otherwise statistics are written here. may be null
   * @param reference must be non-null for CRAM support
   * @param inHeader use this header instead of any found in the input SAM/BAM files (useful if file has no header, or needs changing).
   * @param outHeader use this header on output files (or null to skip writing any header)
   * @param terminateBlockedGzip true if BAM or block compressed SAM should have terminator block
   * @throws java.io.IOException if an IO error occurs
   */
  public void mergeSamFiles(Collection<File> samFiles, Collection<File> calibrationFiles, File output, OutputStream out, SequencesReader reference, SAMFileHeader inHeader, SAMFileHeader outHeader, boolean terminateBlockedGzip) throws IOException {
    final boolean isStdio = FileUtils.isStdio(output);
    if (!isStdio) {
      if (calibrationFiles.size() > 0) {
        if (calibrationFiles.size() != samFiles.size()) {
          Diagnostic.warning("Number of calibration files does not match number of SAM files, will not merge calibration files.");
        } else if (mFilterParams.isFiltering()) {
          Diagnostic.warning("Depending on filter criteria, output calibration data may be inaccurate. Consider recalibrating output file.");
        }
      }
    }
    final File alignmentOutputFile;
    final long recordsIn;
    final long recordsOut;
    final SingletonPopulatorFactory<SamReaderRecord> pf = new SingletonPopulatorFactory<>(new SamReaderRecord.SamReaderRecordPopulator());
    final SamReadingContext context = new SamReadingContext(samFiles, mNumberThreads, mFilterParams, inHeader, reference);
    try (final ThreadedMultifileIterator<SamReaderRecord> outer = new ThreadedMultifileIterator<>(context, pf)) {
      final RecordIterator<SamReaderRecord> it = context.filterParams().findAndRemoveDuplicates() ? new DedupifyingRecordIterator<>(outer) : outer;
      if (outHeader != null) {
        outHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        if (mAddProgramRecord) {
          SamUtils.addProgramRecord(outHeader);
        }
        SamUtils.updateRunId(outHeader);
      }

      try (SamOutput so = SamOutput.getSamOutput(output, out, outHeader == null ? inHeader : outHeader, mGzip, true, outHeader != null, terminateBlockedGzip, mCreateIndex, reference)) {
        alignmentOutputFile = so.getOutFile();
        try (SAMFileWriter writer = so.getWriter()) {
          while (it.hasNext()) {
            final SAMRecord rec = it.next().inner();
            if (mLegacy) {
              SamUtils.convertToLegacyCigar(rec);
            }
            if (mRenameWithRg) {
              renameWithReadGroupId(rec);
            }
            writer.addAlignment(rec);
          }
        }
      }
      recordsIn = it.getTotalRecordsCount();
      recordsOut = it.getTotalRecordsCount() - it.getFilteredRecordsCount() - it.getDuplicateRecordsCount() - it.getInvalidRecordsCount();
    }
    if (!isStdio) {
      if (calibrationFiles.size() > 0 && calibrationFiles.size() == samFiles.size()) {
        final Calibrator c = new Calibrator(Calibrator.getCovariateSet(calibrationFiles.iterator().next()), null);
        for (final File f : calibrationFiles) {
          c.accumulate(f);
        }
        c.writeToFile(new File(alignmentOutputFile.getParent(), alignmentOutputFile.getName() + CommonFlags.RECALIBRATE_EXTENSION));
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

  private static void renameWithReadGroupId(SAMRecord rec) {
    final SAMReadGroupRecord readGroup = rec.getReadGroup();
    if (readGroup == null || readGroup.getId() == null) {
      return;
    }
    rec.setReadName(readGroup.getId() + "-" + rec.getReadName());
  }
}
