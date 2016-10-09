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
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import com.rtg.bed.BedReader;
import com.rtg.bed.BedRecord;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.sam.RecordIterator;
import com.rtg.sam.SamOutput;
import com.rtg.sam.SamRecordPopulator;
import com.rtg.sam.SamUtils;
import com.rtg.sam.SmartSamWriter;
import com.rtg.sam.ThreadedMultifileIterator;
import com.rtg.util.MathUtils;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.TextTable;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogStream;

import htsjdk.samtools.SAMRecord;

/**
 * Remove probes from mapped reads using a bed file to describe probes.
 * <ul>
 *   <li>Assumes probe is on first read in sequencing</li>
 *   <li>Assumes 6th column is +/- representing whether probe is on forward or reverse strand</li>
 *   <li>Only updates records that are mapped, and marked as first in sequencing</li>
 *   <li>Update is made by removing values as appropriate from: the read string, the quality string, and the cigar. As well as updating the alignment start position if required.</li>
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
  private long mTotalMappedReads;
  private long mTotalMappedPos;
  private long mTotalMappedNeg;
  private long mTotalStrippedPos;
  private long mTotalStrippedNeg;
  private long mTotalStrippedReads;


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
    final PosChecker posChecker = new PosChecker(tolerance);
    final NegChecker negChecker = new NegChecker(tolerance);
    final File bedFile = (File) mFlags.getValue(PROBE_BED);

    mTotalMappedReads = 0;
    mTotalMappedPos = 0;
    mTotalMappedNeg = 0;
    mTotalStrippedPos = 0;
    mTotalStrippedNeg = 0;
    mTotalStrippedReads = 0;

    long totalUnmappedReads = 0;
    long totalRecords = 0;
    long totalAmbiguousRecords = 0;
    final HashMap<String, ReadStatus> ambiguousReads = new HashMap<>();

    try (RecordIterator<SAMRecord> reader = new ThreadedMultifileIterator<>(Collections.singletonList(input), new SingletonPopulatorFactory<>(new SamRecordPopulator()))) {
      loadProbeBed(bedFile);
      final File outputFile = new File(outputDirectory(), ALIGNMENT_FILE_NAME);
      try (SamOutput samOutput = SamOutput.getSamOutput(outputFile, out, reader.header(), !mFlags.isSet(CommonFlags.NO_GZIP), true, null)) {
        try (SmartSamWriter writer = new SmartSamWriter(samOutput.getWriter())) {
          while (reader.hasNext()) {
            final SAMRecord record = reader.next();
            if (record.getFirstOfPairFlag()) {
              totalRecords++;
              final int ih = MathUtils.unboxNatural(SamUtils.getNHOrIH(record));
              final boolean unique = 1 == ih || !record.getNotPrimaryAlignmentFlag();
              boolean stripped = false;
              boolean mapped = false;
              if (!record.getReadUnmappedFlag()) {
                mapped = true;

                stripped = checkList(mPosRanges.get(record.getReferenceName()), record, null, tolerance, posChecker)
                  || checkList(mNegRanges.get(record.getReferenceName()), record, null, tolerance, negChecker);
                if (!stripped) {
                  record.setAttribute("XS", "failed");
                }
              }

              final boolean negative = record.getReadNegativeStrandFlag();
              if (unique) {
                if (mapped) {
                  addMapped(stripped, negative);
                } else {
                  totalUnmappedReads++;
                }
              } else {
                // If not uniquely mapped, accumulate and delay addition to the stats until the end
                totalAmbiguousRecords++;
                ReadStatus s = ambiguousReads.get(record.getReadName());
                if (s == null) {
                  s = new ReadStatus(record.getReadName(), stripped, negative);
                  ambiguousReads.put(record.getReadName(), s);
                } else if (stripped && !s.stripped()) {
                  s.setStatus(stripped, negative);
                }
              }
            }
            writer.addRecord(record);
          }
          Diagnostic.progress("Closing SAM writer");
        }
      }
    }

    final long totalReads = mTotalMappedReads + totalUnmappedReads + ambiguousReads.size();
    Diagnostic.userLog("Total R1 records: " + totalRecords);
    Diagnostic.userLog("Total R1 reads: " + totalReads);
    Diagnostic.userLog("Total R1 unmapped reads (records): " + totalUnmappedReads);
    Diagnostic.userLog("Total R1 uniquely mapped reads (records): " + mTotalMappedReads);
    Diagnostic.userLog("Total R1 ambiguously mapped records: " + totalAmbiguousRecords);
    Diagnostic.userLog("Total R1 ambiguously mapped reads: " + ambiguousReads.size());

    // Now merge in counts of any ambiguously mapped reads
    for (ReadStatus s : ambiguousReads.values()) {
      addMapped(s.stripped(), s.negative());
    }

    Diagnostic.userLog("PROBE OFFSETS");
    final TextTable offsetSummary = new TextTable();
    offsetSummary.addRow("Delta", "+", "-");
    offsetSummary.addSeparator();
    for (int i = 0; i < posChecker.mPosDiffStats.length; i++) {
      offsetSummary.addRow(Integer.toString(i - tolerance), Integer.toString(posChecker.mPosDiffStats[i]), Integer.toString(negChecker.mPosDiffStats[tolerance * 2 - i]));
    }
    Diagnostic.userLog(offsetSummary.toString());
    FileUtils.stringToFile(offsetSummary.getAsTsv(), new File(outputDirectory(), PROBE_OFFSET_TABLE_FILE));

    Diagnostic.userLog("CIGAR OPERATIONS WITHIN PROBE");
    final TextTable cigarSummary = new TextTable();
    cigarSummary.addRow("Strand", "Length", "X", "I", "D", "S");
    cigarSummary.addSeparator();
    for (int i = 0; i < PositionAndStrandChecker.MAX_OP_LEN; i++) {
      addCigarRow(cigarSummary, "+", posChecker, i);
      addCigarRow(cigarSummary, "-", negChecker, i);
    }
    Diagnostic.userLog(cigarSummary.toString());
    FileUtils.stringToFile(cigarSummary.getAsTsv(), new File(outputDirectory(), CIGAR_OP_TABLE_FILE));

    final long totalstripped = mTotalStrippedPos + mTotalStrippedNeg;
    final long total = mTotalMappedPos + mTotalMappedNeg;
    final TextTable summary = new TextTable();
    summary.addRow("Strand", "Alignments", "Stripped", "Identified", "nt/read");
    summary.addSeparator();
    summary.addRow("+", Long.toString(mTotalMappedPos), Long.toString(mTotalStrippedPos), String.format("%.2f%%", (double) mTotalStrippedPos / mTotalMappedPos * 100.0), String.format("%.1f", (double) posChecker.mBasesTrimmed / mTotalStrippedPos));
    summary.addRow("-", Long.toString(mTotalMappedNeg), Long.toString(mTotalStrippedNeg), String.format("%.2f%%", (double) mTotalStrippedNeg / mTotalMappedNeg * 100.0), String.format("%.1f", (double) negChecker.mBasesTrimmed / mTotalStrippedNeg));
    summary.addRow("Both", Long.toString(total), Long.toString(totalstripped), String.format("%.2f%%", (double) totalstripped / total * 100.0), String.format("%.1f", (double) (posChecker.mBasesTrimmed + negChecker.mBasesTrimmed) / totalstripped));
    FileUtils.stringToFile(summary.getAsTsv(), new File(outputDirectory(), PROBE_SUMMARY_FILE));

    Diagnostic.userLog("ON TARGET RATES");
    final TextTable onTargetSummary = new TextTable();
    onTargetSummary.addRow("Group", "Total", "Mapped", "On Target", "%/Total", "%/Mapped");
    onTargetSummary.addRow("Reads", Long.toString(totalReads), Long.toString(mTotalMappedReads), Long.toString(mTotalStrippedReads), String.format("%.2f%%", (double) mTotalStrippedReads / (double) totalReads * 100.0), String.format("%.2f%%", (double) mTotalStrippedReads / (double) mTotalMappedReads * 100.0));
    Diagnostic.userLog(onTargetSummary.toString());
    FileUtils.stringToFile(onTargetSummary.getAsTsv(), new File(outputDirectory(), ON_TARGET_SUMMARY_FILE));

    try (PrintStream summaryOut = new PrintStream(FileUtils.createTeedOutputStream(FileUtils.createOutputStream(new File(outputDirectory(), CommonFlags.SUMMARY_FILE), false), out))) {
      summaryOut.println(onTargetSummary);
    }
    return 0;
  }

  private void addMapped(boolean stripped, boolean negative) {
    mTotalMappedReads++;
    if (stripped) {
      mTotalStrippedReads++;
    }
    if (negative) {
      mTotalMappedNeg++;
      if (stripped) {
        mTotalStrippedNeg++;
      }
    } else {
      mTotalMappedPos++;
      if (stripped) {
        mTotalStrippedPos++;
      }
    }
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
    } else {
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

  private static final class ReadStatus {
    final String mName;
    boolean mStripped = false;
    boolean mNegative = false;
    private ReadStatus(String name, boolean stripped, boolean negative) {
      mName = name;
      setStatus(stripped, negative);
    }
    public int hashCode() {
      return mName.hashCode();
    }
    public boolean equals(Object that) {
      return that instanceof ReadStatus && mName.equals(((ReadStatus) that).mName);
    }
    private boolean stripped() {
      return mStripped;
    }
    private boolean negative() {
      return mNegative;
    }
    private void setStatus(boolean stripped, boolean negative) {
      mStripped = stripped;
      mNegative = negative;
    }
  }
}
