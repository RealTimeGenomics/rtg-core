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
package com.rtg.sam.probe;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import com.rtg.bed.BedReader;
import com.rtg.bed.BedRecord;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.sam.SamOutput;
import com.rtg.sam.SortingSamReader;
import com.rtg.util.TextTable;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogStream;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

/**
 * Remove probes from mapped reads using a bed file to describe probes.
 * <ul>
 *   <li>Assumes probe is on first read in sequencing</li>
 *   <li>Assumes 6th column is +/- representing whether probe is on forward or reverse strand</li>
 *   <li>Only updates records that are mapped, and marked as first in sequencing</li>
 *   <li>Update is made by removing values as appropriate from: the read string, the quality string, and the cigar. As well as updating the alignment start position if required.</li>
 *   <li>No changes are made to the mate position or <code>TLEN</code> field of any record (so these are likely to be incorrect when records have been updated)</li>
 *   <li>Mapped first arm reads that do not match a known probe will receive an annotation (XS:Z:failed)</li>
 * </ul>
 */
public class DeProbeCli extends LoggedCli {

  private static final String PROBE_BED = "probe-bed";
  private static final String TOLERANCE_FLAG = "tolerance";
  static final String ALIGNMENT_FILE_NAME = "alignments_deprobed.bam";
  static final String PROBE_OFFSET_TABLE_FILE = "probe_offsets.tsv";
  static final String CIGAR_OP_TABLE_FILE = "cigar_ops.tsv";
  static final String PROBE_SUMMARY_FILE = "probe_summary.tsv";
  static final String ON_TARGET_SUMMARY_FILE = "on_target.tsv";


  private ReferenceRanges<String> mPosRanges;
  private ReferenceRanges<String> mNegRanges;


  @Override
  public String moduleName() {
    return "deprobe";
  }

  @Override
  public String description() {
    return "trim strand-specific probes from mapped SAM/BAM alignments";
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  private static void initFlags(CFlags flags) {
    flags.registerExtendedHelp();
    CommonFlagCategories.setCategories(flags);
    flags.registerRequired('i', CommonFlags.INPUT_FLAG, File.class, "FILE", "SAM/BAM file containing input alignments").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerRequired('o', CommonFlags.OUTPUT_FLAG, File.class, "FILE", "directory for output").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerRequired('b', PROBE_BED, File.class, "FILE", "BED file specifying each probe location and strand").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerOptional(TOLERANCE_FLAG, Integer.class, CommonFlags.INT, "start position tolerance for probe matching", 3).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  @Override
  protected int mainExec(OutputStream out, LogStream log) throws IOException {
    final int tolerance = (Integer) mFlags.getValue(TOLERANCE_FLAG);
    final File input = (File) mFlags.getValue(CommonFlags.INPUT_FLAG);
    long totalPos = 0;
    long totalNeg = 0;
    long posStripped = 0;
    long negStripped = 0;
    long totalReads = 0;
    long onTargetReads = 0;
    long totalRecords = 0;
    long onTargetRecords = 0;
    long totalMappedBases = 0;
    long totalMappedRecords = 0;
    long totalMappedReads = 0;
    long totalBases = 0;
    long onTargetBases = 0;
    final PosChecker posChecker = new PosChecker(tolerance);
    final NegChecker negChecker = new NegChecker(tolerance);
    final File bedFile = (File) mFlags.getValue(PROBE_BED);
    final File outputDir = (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);

    Diagnostic.progress("Beginning read name sort");
    try (SortingSamReader reader = SortingSamReader.reader(input, SAMFileHeader.SortOrder.queryname, outputDir)) {
      Diagnostic.progress("read name sort complete");
      loadProbeBed(bedFile);
//      final ReferenceRanges<String> bedReferenceRanges = SamRangeUtils.createBedReferenceRanges(reader.getFileHeader(), (File) mFlags.getValue(PROBE_BED));
      final File outputFile = new File(outputDir, ALIGNMENT_FILE_NAME);
      final SAMFileHeader fileHeader = reader.getFileHeader();
      fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
      try (SamOutput samOutput = SamOutput.getSamOutput(outputFile, out, fileHeader, !mFlags.isSet(CommonFlags.NO_GZIP), false)) {
        try (SAMFileWriter writer = samOutput.getWriter()) {
          List<SAMRecord> readRecords;
          final ReadRecordGrouper grouper = new ReadRecordGrouper(reader);
          final Map<Integer, SAMRecord> posMap = new HashMap<>();
          while ((readRecords = grouper.nextRead()) != null) {
            posMap.clear();
            totalReads++;
            for (SAMRecord record : readRecords) {
              totalRecords++;
              totalBases += record.getReadLength();
              posMap.put(record.getAlignmentStart(), record);
            }
            boolean strippedRead = false;
            boolean mappedRead = false;
            for (SAMRecord record : readRecords) {
              final int initialLength = record.getReadLength();
              if (!record.getReadUnmappedFlag()) {
                mappedRead = true;
                totalMappedRecords++;
                totalMappedBases += record.getReadLength();
                final SAMRecord mate = posMap.get(record.getMateAlignmentStart());
                if (record.getFirstOfPairFlag()) {
                  if (checkList(mPosRanges.get(record.getReferenceName()), record, mate, tolerance, posChecker)
                    || checkList(mNegRanges.get(record.getReferenceName()), record, mate, tolerance, negChecker)) {
                    strippedRead = true;
                    onTargetRecords++;
                    onTargetBases += initialLength;
                    if (record.getProperPairFlag() && mate != null) {
                      onTargetRecords++;
                      onTargetBases += mate.getReadLength();
                    }
                    if (record.getReadNegativeStrandFlag()) {
                      negStripped++;
                    } else {
                      posStripped++;
                    }
                  } else {
                    record.setAttribute("XS", "failed");
                  }
                  if (record.getReadNegativeStrandFlag()) {
                    totalNeg++;
                  } else {
                    totalPos++;
                  }
                }
              }
              writer.addAlignment(record);
            }
            if (strippedRead) {
              onTargetReads++;
            }
            if (mappedRead) {
              totalMappedReads++;
            }
          }
          Diagnostic.progress("Closing SAM writer");
        }
      }
    }

    Diagnostic.userLog("PROBE OFFSETS");
    final TextTable offsetSummary = new TextTable();
    offsetSummary.addRow("Delta", "+", "-");
    offsetSummary.addSeparator();
    for (int i = 0; i < posChecker.mPosDiffStats.length; i++) {
      offsetSummary.addRow(Integer.toString(i - tolerance), Integer.toString(posChecker.mPosDiffStats[i]), Integer.toString(negChecker.mPosDiffStats[tolerance * 2 - i]));
    }
    Diagnostic.userLog(offsetSummary.toString());
    FileUtils.stringToFile(offsetSummary.getAsTsv(), new File(outputDir, PROBE_OFFSET_TABLE_FILE));

    Diagnostic.userLog("CIGAR OPERATIONS WITHIN PROBE");
    final TextTable cigarSummary = new TextTable();
    cigarSummary.addRow("Strand", "Length", "X", "I", "D", "S");
    cigarSummary.addSeparator();
    for (int i = 0; i < PositionAndStrandChecker.MAX_OP_LEN; i++) {
      addCigarRow(cigarSummary, "+", posChecker, i);
      addCigarRow(cigarSummary, "-", negChecker, i);
    }
    Diagnostic.userLog(cigarSummary.toString());
    FileUtils.stringToFile(cigarSummary.getAsTsv(), new File(outputDir, CIGAR_OP_TABLE_FILE));

    final TextTable onTargetSummary = new TextTable();
    onTargetSummary.addRow("Group", "Total", "Mapped", "On Target", "%/Total", "%/Mapped");
    onTargetSummary.addRow("Reads", Long.toString(totalReads), Long.toString(totalMappedReads), Long.toString(onTargetReads), String.format("%.2f%%", (double) onTargetReads / (double) totalReads * 100.0), String.format("%.2f%%", (double) onTargetReads / (double) totalMappedReads * 100.0));
    onTargetSummary.addRow("Records", Long.toString(totalRecords), Long.toString(totalMappedRecords), Long.toString(onTargetRecords), String.format("%.2f%%", (double) onTargetRecords / (double) totalRecords * 100.0), String.format("%.2f%%", (double) onTargetRecords / (double) totalMappedRecords * 100.0));
    onTargetSummary.addRow("Bases", Long.toString(totalBases), Long.toString(totalMappedBases), Long.toString(onTargetBases), String.format("%.2f%%", (double) onTargetBases / (double) totalBases * 100.0), String.format("%.2f%%", (double) onTargetBases / (double) totalMappedBases * 100.0));
    Diagnostic.userLog(onTargetSummary.toString());
    FileUtils.stringToFile(onTargetSummary.getAsTsv(), new File(outputDir, ON_TARGET_SUMMARY_FILE));

    try (PrintStream summaryOut = new PrintStream(FileUtils.createTeedOutputStream(FileUtils.createOutputStream(new File(outputDir, CommonFlags.SUMMARY_FILE), false), out))) {
      final long totalstripped = posStripped + negStripped;
      final long total = totalPos + totalNeg;
      final TextTable summary = new TextTable();
      summary.addRow("Strand", "Alignments", "Stripped", "Identified", "nt/read");
      summary.addSeparator();
      summary.addRow("+", Long.toString(totalPos), Long.toString(posStripped), String.format("%.2f%%", (double) posStripped / totalPos * 100.0), String.format("%.1f", (double) posChecker.mBasesTrimmed / posStripped));
      summary.addRow("-", Long.toString(totalNeg), Long.toString(negStripped), String.format("%.2f%%", (double) negStripped / totalNeg * 100.0), String.format("%.1f", (double) negChecker.mBasesTrimmed / negStripped));
      summary.addRow("Both", Long.toString(total), Long.toString(totalstripped), String.format("%.2f%%", (double) totalstripped / total * 100.0), String.format("%.1f", (double) (posChecker.mBasesTrimmed + negChecker.mBasesTrimmed) / totalstripped));
      summaryOut.println(summary);
      FileUtils.stringToFile(summary.getAsTsv(), new File(outputDir, PROBE_SUMMARY_FILE));
    }
    return 0;
  }


  private void addCigarRow(TextTable ctype, String strand, PositionAndStrandChecker checker, int i) {
    final String len = Integer.toString(i + 1) + (i == PositionAndStrandChecker.MAX_OP_LEN - 1 ? "+" : "");
    ctype.addRow(strand, len, Integer.toString(checker.mMismatchStats[i]), Integer.toString(checker.mInsertStats[i]), Integer.toString(checker.mDeletionStats[i]), Integer.toString(checker.mSoftClipStats[i]));
  }

  private void loadProbeBed(File bedFile) throws IOException {
    final ReferenceRanges.Accumulator<String> posRangesAccum = new ReferenceRanges.Accumulator<>();
    final ReferenceRanges.Accumulator<String> negRangesAccum = new ReferenceRanges.Accumulator<>();
    try (BedReader bedReader = BedReader.openBedReader(null, bedFile, 3)) {
      while (bedReader.hasNext()) {
        final BedRecord record = bedReader.next();
        final RangeList.RangeData<String> rangeData = new RangeList.RangeData<>(record.getStart(), record.getEnd(), Arrays.asList(record.getAnnotations()));
        if ("+".equals(record.getAnnotations()[2])) {
          posRangesAccum.addRangeData(record.getSequenceName(), rangeData);
        } else if ("-".equals(record.getAnnotations()[2])) {
          negRangesAccum.addRangeData(record.getSequenceName(), rangeData);
        } else {
          throw new IOException("BED record does not indicate strand as '-' or '+' in column 6: " + record);
        }
      }
    }
    mPosRanges = posRangesAccum.getReferenceRanges();
    mNegRanges = negRangesAccum.getReferenceRanges();
  }

  private static boolean checkList(RangeList<String> list, SAMRecord record, SAMRecord mate, int tolerance, PositionAndStrandChecker checker) {
    if (list == null) {
      return false;
    }
    int index = checker.getStartDataIndex(record, list);
    RangeList.RangeData<String> data = list.getFullRangeList().get(index);
    while (recordOverlap(record, data, tolerance)) {
      final List<String> metas = data.getMeta();
      if (metas != null && metas.size() >= 3) {
        if (checker.check(record, data)) {
          checker.stripRecord(record, mate, data);
          return true;
        }
      }
      if (list.getFullRangeList().size() > ++index) {
        data = list.getFullRangeList().get(index);
      } else {
        data = null;
      }
    }
    return false;
  }

  private static boolean recordOverlap(SAMRecord rec, RangeList.RangeData<String> data, int tolerance) {
    if (data == null) {
      return false;
    }
    final int dataMin;
    if (data.getStart() - tolerance > data.getStart()) {
      dataMin = Integer.MIN_VALUE;
    } else {
      dataMin = data.getStart() - tolerance;
    }
    final int dataMax;
    if (data.getEnd() + tolerance < data.getEnd()) {
      dataMax = Integer.MAX_VALUE;
    } else  {
      dataMax = data.getEnd() + tolerance;
    }
    final int alignmentStart = rec.getAlignmentStart() - 1;
    if (alignmentStart >= dataMin && alignmentStart < dataMax) {
      return true;
    }
    final int alignmentEnd = rec.getAlignmentEnd();
    if (alignmentEnd > dataMin && alignmentEnd < dataMax) {
      return true;
    }
    return false;
  }

  private static class ReadRecordGrouper {
    private final Iterator<SAMRecord> mReader;
    private SAMRecord mCurrent = null;

    ReadRecordGrouper(Iterable<SAMRecord> reader) {
      mReader = reader.iterator();
    }

    void initCurrent() {
      if (mCurrent == null && mReader.hasNext()) {
        mCurrent = mReader.next();
      }
    }

    List<SAMRecord> nextRead() {
      initCurrent();
      if (mCurrent == null) {
        return null;
      }
      final ArrayList<SAMRecord> read = new ArrayList<>();
      read.add(mCurrent);
      final String name = mCurrent.getReadName();
      while (true) {
        if (mReader.hasNext()) {
          mCurrent = mReader.next();
          if (mCurrent.getReadName().equals(name)) {
            read.add(mCurrent);
          } else {
            break;
          }
        } else {
          mCurrent = null;
          break;
        }
      }
      return read;
    }
  }

}
