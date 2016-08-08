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
import java.util.List;

import com.rtg.bed.BedReader;
import com.rtg.bed.BedRecord;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.sam.SamOutput;
import com.rtg.sam.SamUtils;
import com.rtg.sam.SmartSamWriter;
import com.rtg.util.TextTable;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.ReferenceRanges;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

/**
 *
 */
public class BamStripProbes extends AbstractCli {

  private static final String PROBE_BED = "probe-bed";
  private static final String TOLERANCE_FLAG = "tolerance";


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
    flags.registerRequired('o', CommonFlags.OUTPUT_FLAG, File.class, "FILE", "name for output SAM/BAM file").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerRequired('b', PROBE_BED, File.class, "FILE", "BED file specifying each probe location and strand").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerOptional(TOLERANCE_FLAG, Integer.class, CommonFlags.INT, "start position tolerance for probe matching", 3).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    CommonFlags.initNoGzip(flags);
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final int tolerance = (Integer) mFlags.getValue(TOLERANCE_FLAG);
    final File input = (File) mFlags.getValue(CommonFlags.INPUT_FLAG);
    long totalPos = 0;
    long totalNeg = 0;
    long posStripped = 0;
    long negStripped = 0;
    final PosChecker posChecker = new PosChecker(tolerance);
    final NegChecker negChecker = new NegChecker(tolerance);
    final File bedFile = (File) mFlags.getValue(PROBE_BED);
    try (SamReader reader = SamUtils.makeSamReader(input)) {
      loadProbeBed(bedFile);
//      final ReferenceRanges<String> bedReferenceRanges = SamRangeUtils.createBedReferenceRanges(reader.getFileHeader(), (File) mFlags.getValue(PROBE_BED));
      try (SamOutput samOutput = SamOutput.getSamOutput((File) mFlags.getValue(CommonFlags.OUTPUT_FLAG), out, reader.getFileHeader(), !mFlags.isSet(CommonFlags.NO_GZIP))) {
        try (SmartSamWriter writer = new SmartSamWriter(samOutput.getWriter())) {
          for (SAMRecord record : reader) {
            if (!record.getReadUnmappedFlag()) {
              if (record.getFirstOfPairFlag()) {
                if (checkList(mPosRanges.get(record.getReferenceName()), record, tolerance, posChecker)) {
                  posStripped++;
                } else if (checkList(mNegRanges.get(record.getReferenceName()), record, tolerance, negChecker)) {
                  negStripped++;
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
            writer.addRecord(record);
          }
        }
      }
    }

    err.println("PROBE OFFSETS");
    final TextTable offsetSummary = new TextTable();
    offsetSummary.addRow("Delta", "+", "-");
    offsetSummary.addSeparator();
    for (int i = 0; i < posChecker.mPosDiffStats.length; i++) {
      offsetSummary.addRow(Integer.toString(i - tolerance), Integer.toString(posChecker.mPosDiffStats[i]), Integer.toString(negChecker.mPosDiffStats[tolerance * 2 - i]));
    }
    err.println(offsetSummary);

    err.println("CIGAR OPERATIONS WITHIN PROBE");
    final TextTable cigarSummary = new TextTable();
    cigarSummary.addRow("Strand", "Length", "X", "I", "D", "S");
    cigarSummary.addSeparator();
    for (int i = 0; i < PositionAndStrandChecker.MAX_OP_LEN; i++) {
      addCigarRow(cigarSummary, "+", posChecker, i);
      addCigarRow(cigarSummary, "-", negChecker, i);
    }
    err.println(cigarSummary);

    err.println("SUMMARY");
    final long totalstripped = posStripped + negStripped;
    final long total = totalPos + totalNeg;
    final TextTable summary = new TextTable();
    summary.addRow("Strand", "Alignments", "Stripped", "Identified", "nt/read");
    summary.addSeparator();
    summary.addRow("+", Long.toString(totalPos), Long.toString(posStripped), String.format("%.2f%%", (double) posStripped / totalPos * 100.0), String.format("%.1f", (double) posChecker.mBasesTrimmed / posStripped));
    summary.addRow("-", Long.toString(totalNeg), Long.toString(negStripped), String.format("%.2f%%", (double) negStripped / totalNeg * 100.0), String.format("%.1f", (double) negChecker.mBasesTrimmed / posStripped));
    summary.addRow("Both", Long.toString(total), Long.toString(totalstripped), String.format("%.2f%%", (double) totalstripped / total * 100.0), String.format("%.1f", (double) (posChecker.mBasesTrimmed + negChecker.mBasesTrimmed) / totalstripped));
    err.println(summary);

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

  private static boolean checkList(RangeList<String> list, SAMRecord record, int tolerance, PositionAndStrandChecker checker) {
    int index = checker.getStartDataIndex(record, list);
    RangeList.RangeData<String> data = list.getFullRangeList().get(index);
    while (recordOverlap(record, data, tolerance)) {
      final List<String> metas = data.getMeta();
      if (metas != null && metas.size() >= 3) {
        if (checker.check(record, data)) {
          checker.stripRecord(record, data);
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

}
