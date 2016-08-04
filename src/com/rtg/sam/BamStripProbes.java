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
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.rtg.bed.BedReader;
import com.rtg.bed.BedRecord;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.io.FileUtils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
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
    final SamBamBaseFile baseFile = SamBamBaseFile.getBaseFile((File) mFlags.getValue(CommonFlags.OUTPUT_FLAG), !mFlags.isSet(CommonFlags.NO_GZIP));
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
      try (SmartSamWriter writer = new SmartSamWriter(getSAMFileWriter(baseFile, reader))) {
        for (SAMRecord record : reader) {
          if (!record.getReadUnmappedFlag()) {
            if (record.getFirstOfPairFlag()) {
              boolean stripped = false;
              boolean neg = false;
              if (checkList(posRanges.get(record.getReferenceName()), record, tolerance, posChecker)) {
                stripped = true;
              } else if (checkList(negRanges.get(record.getReferenceName()), record, tolerance, negChecker)) {
                stripped = true;
                neg = true;
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
    System.err.println("total = " + total);
    System.err.println("totalstripped = " + totalstripped);
    System.err.println(String.format("percentage Stripped: %.2f%%", (double) totalstripped / total * 100.0));
    System.err.println("totalForward = " + totalForward);
    System.err.println("totalReverse = " + totalReverse);
    System.err.println("forward = " + forward);
    System.err.println("reverse = " + reverse);
    System.err.print(String.format("%9s | %9s %9s%n", "DIFF", "POS", "NEG"));
    for (int i = 0; i < posChecker.mStats.length; i++) {
      System.err.print(String.format("%9d | %9d %9d%n", i - tolerance, posChecker.mStats[i], negChecker.mStats[tolerance * 2 - i]));
    }

    return 0;
  }

  abstract static class PositionAndStrandChecker {
    protected final int mTolerance;
    protected final int[] mStats;
    PositionAndStrandChecker(int tolerance) {
      mTolerance = tolerance;
      mStats = new int[mTolerance * 2 + 1];
    }

    abstract boolean check(SAMRecord record, RangeList.RangeData<String> data);
    abstract int getStartDataIndex(SAMRecord record, RangeList<String> list);
    abstract void stripRecord(SAMRecord record, RangeList.RangeData<String> data);
  }

  private static class PosChecker extends PositionAndStrandChecker {
    PosChecker(int tolerance) {
      super(tolerance);
    }

    @Override
    public boolean check(SAMRecord record, RangeList.RangeData<String> data) {
      final int alignmentStart = record.getAlignmentStart() - 1;
      if (!record.getReadNegativeStrandFlag()) {
        if (data.getStart() > alignmentStart - mTolerance && data.getStart() < alignmentStart + mTolerance) {
//                    System.err.println(record.getSAMString() + " strip forward to: " + data.getEnd() + " (" + data.getStart() + " : " + data.getEnd() + ")");
          return true;
        }
      }
      return false;
    }

    @Override
    public int getStartDataIndex(SAMRecord record, RangeList<String> list) {
      final int alignmentStart = record.getAlignmentStart() - 1;
      return list.findFullRangeIndex(alignmentStart - mTolerance);
    }

    @Override
    void stripRecord(SAMRecord record, RangeList.RangeData<String> data) {
      final int diff = record.getAlignmentStart() - 1 - data.getStart();
      mStats[mTolerance + diff]++;
      setAlignmentStart(record, data.getStart());
    }
  }
  private static class NegChecker extends PositionAndStrandChecker {

    NegChecker(int tolerance) {
      super(tolerance);
    }

    @Override
    public boolean check(SAMRecord record, RangeList.RangeData<String> data) {
      if (record.getReadNegativeStrandFlag()) {
        final int alignmentEnd = record.getAlignmentEnd();
        if (data.getEnd() > alignmentEnd - mTolerance && data.getEnd() < alignmentEnd + mTolerance) {
//                    System.err.println(record.getSAMString() + "strip back to: " + data.getStart() + " (" + data.getStart() + " : " + data.getEnd() + ")");
          return true;
        }
      }
      return false;
    }

    @Override
    public int getStartDataIndex(SAMRecord record, RangeList<String> list) {
      return list.findFullRangeIndex(record.getAlignmentEnd() - mTolerance);
    }

    @Override
    void stripRecord(SAMRecord record, RangeList.RangeData<String> data) {
      final int diff = record.getAlignmentEnd() - data.getEnd();
      mStats[mTolerance + diff]++;
      setAlignmentEnd(record, data.getEnd());
    }
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
    final int alignmentStart = rec.getAlignmentStart() - 1;
    if (alignmentStart >= data.getStart() - tolerance && alignmentStart < data.getEnd() + tolerance) {
      return true;
    }
    final int alignmentEnd = rec.getAlignmentEnd();
    if (alignmentEnd > data.getStart() - tolerance && alignmentEnd < data.getEnd() + tolerance) {
      return true;
    }
    return false;
  }

  private SAMFileWriter getSAMFileWriter(SamBamBaseFile baseFile, SamReader reader) throws IOException {
    if (baseFile.isBam()) {
      return new SAMFileWriterFactory().makeBAMWriter(reader.getFileHeader(), true, FileUtils.createOutputStream(baseFile.suffixedFile(""), false));
    } else {
      return new SAMFileWriterFactory().makeSAMWriter(reader.getFileHeader(), true, FileUtils.createOutputStream(baseFile.suffixedFile(""), baseFile.isGzip()));
    }
  }

  static void setAlignmentStart(SAMRecord record, int alignmentStart) {
    int readStart = 0;
    int refStart = record.getAlignmentStart() - 1;
    final List<CigarElement> cigarElements = new ArrayList<>();
    for (CigarElement e : record.getCigar().getCigarElements()) {
      final CigarOperator operator = e.getOperator();
      if (alignmentStart > refStart) {
        final int consume = operator.consumesReferenceBases() ? Math.min(alignmentStart - refStart, e.getLength()) : e.getLength();
        if (operator.consumesReferenceBases()) {
          refStart += consume;
        }
        if (operator.consumesReadBases()) {
          readStart += consume;
        }
        if (e.getLength() - consume > 0) {
          cigarElements.add(new CigarElement(e.getLength() - consume, operator));
        }
      } else {
        cigarElements.add(e);
      }
    }
    final byte[] readBases = record.getReadBases();
    record.setReadBases(Arrays.copyOfRange(readBases, readStart, readBases.length));
    final byte[] baseQualities = record.getBaseQualities();
    record.setBaseQualities(Arrays.copyOfRange(baseQualities, readStart, readBases.length));
    record.setCigar(new Cigar(cigarElements));
    record.setAlignmentStart(alignmentStart + 1);
  }

  static void setAlignmentEnd(SAMRecord record, int alignmentEnd) {
    //end positions are 1 based exclusive
    int readEnd = record.getReadLength();
    int refEnd = record.getAlignmentEnd();
    final List<CigarElement> newCigarElements = new ArrayList<>();
    final List<CigarElement> cigarElements = record.getCigar().getCigarElements();
    for (int i = cigarElements.size() - 1; i >= 0; i--) {
      final CigarElement e = cigarElements.get(i);
      final CigarOperator operator = e.getOperator();
      if (alignmentEnd < refEnd) {
        final int consume = operator.consumesReferenceBases() ? Math.min(refEnd - alignmentEnd, e.getLength()) : e.getLength();
        if (operator.consumesReferenceBases()) {
          refEnd -= consume;
        }
        if (operator.consumesReadBases()) {
          readEnd -= consume;
        }
        if (e.getLength() - consume > 0) {
          newCigarElements.add(new CigarElement(e.getLength() - consume, operator));
        }
      } else {
        newCigarElements.add(e);
      }
    }
    Collections.reverse(newCigarElements);
    final byte[] readBases = record.getReadBases();
    record.setReadBases(Arrays.copyOfRange(readBases, 0, readEnd));
    final byte[] baseQualities = record.getBaseQualities();
    record.setBaseQualities(Arrays.copyOfRange(baseQualities, 0, readEnd));
    record.setCigar(new Cigar(newCigarElements));
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
