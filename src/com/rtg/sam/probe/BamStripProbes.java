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

  private static final String NAME = "sam-strip-probes";

  private static final String DESC = "Remove probes from mapped SAM/BAM";
  private static final String PROBE_BED = "probe-bed";
  private static final String TOLERANCE_FLAG = "tolerance";

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  private static void initFlags(CFlags flags) {
    flags.registerExtendedHelp();
    CommonFlagCategories.setCategories(flags);
    flags.registerRequired('i', CommonFlags.INPUT_FLAG, File.class, "FILE", "SAM or BAM input file").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerRequired('o', CommonFlags.OUTPUT_FLAG, File.class, "FILE", "SAM or BAM output file").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerRequired('b', PROBE_BED, File.class, "FILE", "Bed file for probes").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerOptional(TOLERANCE_FLAG, Integer.class, CommonFlags.INT, "Start position tolerance for probe", 3).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    CommonFlags.initNoGzip(flags);
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final int tolerance = (Integer) mFlags.getValue(TOLERANCE_FLAG);
    final File input = (File) mFlags.getValue(CommonFlags.INPUT_FLAG);
    long totalstripped = 0;
    long total = 0;
    long totalForward = 0;
    long totalReverse = 0;
    long forward = 0;
    long reverse = 0;
    final PosChecker posChecker = new PosChecker(tolerance);
    final NegChecker negChecker = new NegChecker(tolerance);
    try (SamReader reader = SamUtils.makeSamReader(input)) {
      final ReferenceRanges.Accumulator<String> posRangesAccum = new ReferenceRanges.Accumulator<>();
      final ReferenceRanges.Accumulator<String> negRangesAccum = new ReferenceRanges.Accumulator<>();
      try (BedReader bedReader = BedReader.openBedReader(null, (File) mFlags.getValue(PROBE_BED), 3)) {
        while (bedReader.hasNext()) {
          final BedRecord record = bedReader.next();
          final RangeList.RangeData<String> rangeData = new RangeList.RangeData<>(record.getStart(), record.getEnd(), Arrays.asList(record.getAnnotations()));
          if ("+".equals(record.getAnnotations()[2])) {
            posRangesAccum.addRangeData(record.getSequenceName(), rangeData);
          } else if ("-".equals(record.getAnnotations()[2])) {
            negRangesAccum.addRangeData(record.getSequenceName(), rangeData);
          } else {
            throw new RuntimeException("badness: " + record);
          }
        }
      }
      final ReferenceRanges<String> posRanges = posRangesAccum.getReferenceRanges();
      final ReferenceRanges<String> negRanges = negRangesAccum.getReferenceRanges();
//      final ReferenceRanges<String> bedReferenceRanges = SamRangeUtils.createBedReferenceRanges(reader.getFileHeader(), (File) mFlags.getValue(PROBE_BED));
      try (SamOutput samOutput = SamOutput.getSamOutput((File) mFlags.getValue(CommonFlags.OUTPUT_FLAG), out, reader.getFileHeader(), !mFlags.isSet(CommonFlags.NO_GZIP))) {
        try (SmartSamWriter writer = new SmartSamWriter(samOutput.getWriter())) {
          for (SAMRecord record : reader) {
            if (!record.getReadUnmappedFlag()) {
              if (record.getFirstOfPairFlag()) {
                boolean stripped = false;
                if (checkList(posRanges.get(record.getReferenceName()), record, tolerance, posChecker)) {
                  stripped = true;
                } else if (checkList(negRanges.get(record.getReferenceName()), record, tolerance, negChecker)) {
                  stripped = true;
                }
                if (record.getReadNegativeStrandFlag()) {
                  totalReverse++;
                } else {
                  totalForward++;
                }
                if (stripped) {
                  totalstripped++;
                  if (record.getReadNegativeStrandFlag()) {
                    reverse++;
                  } else {
                    forward++;
                  }
                } else {
                  record.setAttribute("XS", "failed");
                }
                total++;
              }
            }
            writer.addRecord(record);
          }
        }
      }
    }
    err.println("total = " + total);
    err.println("totalstripped = " + totalstripped);
    err.println(String.format("percentage Stripped: %.2f%%", (double) totalstripped / total * 100.0));
    err.println("totalForward = " + totalForward);
    err.println("totalReverse = " + totalReverse);
    err.println("forward = " + forward);
    err.println("reverse = " + reverse);
    err.print(String.format("%9s | %9s %9s%n", "DIFF", "POS", "NEG"));
    for (int i = 0; i < posChecker.mPosDiffStats.length; i++) {
      err.print(String.format("%9d | %9d %9d%n", i - tolerance, posChecker.mPosDiffStats[i], negChecker.mPosDiffStats[tolerance * 2 - i]));
    }
    err.print(String.format("(S) %4s  | %9s %9s %9s %9s%n", "LEN", "SOFT", "MISM", "INS", "DEL"));
    for (int i = 0; i < PositionAndStrandChecker.MAX_OP_LEN  - 1; i++) {
      err.print(String.format("(+) %4d  | %9d %9d %9d %9d%n", i + 1, posChecker.mSoftClipStats[i], posChecker.mMismatchStats[i], posChecker.mInsertStats[i], posChecker.mDeletionStats[i]));
      err.print(String.format("(-) %4d  | %9d %9d %9d %9d%n", i + 1, negChecker.mSoftClipStats[i], negChecker.mMismatchStats[i], negChecker.mInsertStats[i], negChecker.mDeletionStats[i]));
    }
    final int maxIndex = PositionAndStrandChecker.MAX_OP_LEN - 1;
    err.print(String.format("(+) %4d+ | %9d %9d %9d %9d%n", maxIndex + 1, posChecker.mSoftClipStats[maxIndex], posChecker.mMismatchStats[maxIndex], posChecker.mInsertStats[maxIndex], posChecker.mDeletionStats[maxIndex]));
    err.print(String.format("(-) %4d+ | %9d %9d %9d %9d%n", maxIndex + 1, negChecker.mSoftClipStats[maxIndex], negChecker.mMismatchStats[maxIndex], negChecker.mInsertStats[maxIndex], negChecker.mDeletionStats[maxIndex]));

    return 0;
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

  @Override
  public String description() {
    return DESC;
  }

  @Override
  public String moduleName() {
    return NAME;
  }

}
