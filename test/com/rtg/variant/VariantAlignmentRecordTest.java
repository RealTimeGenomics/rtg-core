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
package com.rtg.variant;

import java.util.Arrays;

import com.rtg.mode.DnaUtils;
import com.rtg.reader.Arm;
import com.rtg.reader.CgUtils;
import com.rtg.reader.FastaUtils;
import com.rtg.sam.ReadGroupUtils;
import com.rtg.sam.SamUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.realign.RealignParams;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import junit.framework.TestCase;

/**
 */
public class VariantAlignmentRecordTest extends TestCase {

  private static final String READ_GROUP_ID = "ASDF";
  private static final String CG1_READ_GROUP_ID = "CG";
  private static final String CG2_READ_GROUP_ID = "CG2";

  public void testBasicOperation() {
    final SAMFileHeader header = new SAMFileHeader();
    final SAMRecord rec = new SAMRecord(header);
    rec.setReadString("TATT");
    rec.setReadPairedFlag(true);
    rec.setAlignmentStart(42);
    rec.setCigarString("2M2D2M");
    rec.setMappingQuality(43);
    rec.setAttribute("NH", 2);
    final byte[] q = {32, 33, 34, 35};
    rec.setBaseQualities(q);
    rec.setInferredInsertSize(-232);
    final VariantAlignmentRecord r = new VariantAlignmentRecord(rec);
    assertEquals("TATT", DnaUtils.bytesToSequenceIncCG(r.getRead()));
    assertEquals(41, r.getStart());
    assertEquals(6, r.getLength());
    assertTrue(Arrays.equals(q, r.getRecalibratedQuality()));
    assertEquals("2M2D2M", r.getCigar());
    assertEquals(43, r.getMappingQuality());
    assertEquals(2, r.getNHOrIH());
    assertEquals("41 2M2D2M TATT ABCD", r.toString());
    assertEquals(false, r.isNegativeStrand());
    assertNull(r.chain());
    r.setNextInChain(r);
    assertEquals(r, r.chain());
    assertFalse(r.isOverflow());
    assertTrue(VariantAlignmentRecord.overflow(0, 1).isOverflow());
    assertFalse(r.isUnmapped());
    assertEquals(-232, r.getFragmentLength());
  }

  public void testNegStrand() {
    final SAMFileHeader header = new SAMFileHeader();
    final SAMRecord rec = new SAMRecord(header);
    rec.setReadNegativeStrandFlag(true);
    VariantAlignmentRecord r = new VariantAlignmentRecord(rec);
    assertTrue(r.isNegativeStrand());
    rec.setReadNegativeStrandFlag(false);
    r = new VariantAlignmentRecord(rec);
    assertFalse(r.isNegativeStrand());
    rec.setReadPairedFlag(true);
    r = new VariantAlignmentRecord(rec);
    assertFalse(r.isNegativeStrand());
    rec.setReadNegativeStrandFlag(true);
    r = new VariantAlignmentRecord(rec);
    assertTrue(r.isNegativeStrand());
    assertFalse(r.isUnmapped());
  }

  public void testUnmapped() {
    final SAMFileHeader header = new SAMFileHeader();
    final SAMRecord rec = new SAMRecord(header);
    rec.setReadUnmappedFlag(true);
    final VariantAlignmentRecord r = new VariantAlignmentRecord(rec);
    assertTrue(r.isUnmapped());
  }

  private void check(final VariantAlignmentRecord a, final VariantAlignmentRecord b) {
    assertTrue(a.compareTo(b) > 0);
    assertTrue(b.compareTo(a) < 0);
    assertEquals(0, a.compareTo(a));
  }

  public void testComparator() {
    final SAMFileHeader header = new SAMFileHeader();
    header.addSequence(new SAMSequenceRecord("seq0", 1000));
    header.addSequence(new SAMSequenceRecord("seq1", 1000));
    header.addSequence(new SAMSequenceRecord("seq2", 1000));
    final SAMRecord rec1 = new SAMRecord(header);
    final SAMRecord rec2 = new SAMRecord(header);
    rec1.setReferenceIndex(-1);
    rec2.setReferenceIndex(1);
    check(new VariantAlignmentRecord(rec1), new VariantAlignmentRecord(rec2));
    rec1.setReferenceIndex(2);
    check(new VariantAlignmentRecord(rec1), new VariantAlignmentRecord(rec2));
    rec1.setReferenceIndex(1);
    rec1.setAlignmentStart(11);
    rec2.setAlignmentStart(10);
    check(new VariantAlignmentRecord(rec1), new VariantAlignmentRecord(rec2));
    rec1.setAlignmentStart(10);
    rec1.setCigarString("11M");
    rec2.setCigarString("10M");
    check(new VariantAlignmentRecord(rec1), new VariantAlignmentRecord(rec2));
    rec2.setCigarString("5M6X");
    check(new VariantAlignmentRecord(rec1), new VariantAlignmentRecord(rec2));
    rec2.setCigarString("11M");
    rec1.setMappingQuality(10);
    rec2.setMappingQuality(11);
    check(new VariantAlignmentRecord(rec1), new VariantAlignmentRecord(rec2));
    rec2.setMappingQuality(10);
    // there are more fields for even finer disambiguation
  }

  static class MachineErrorChooserAdapter implements MachineErrorChooserInterface {
    @Override
    public PhredScaler machineErrors(SAMReadGroupRecord rg, boolean readPaired) {
      return null;
    }

    @Override
    public RealignParams realignParams(SAMReadGroupRecord rg, boolean readPaired) {
      return null;
    }
    @Override
    public MachineType machineType(SAMReadGroupRecord rg, boolean readPaired) {
      return null;
    }

  }

  public void testRecalibrationReadPos() {
    final byte[] q = {32, 33, 34, 35};
    final SAMRecord rec = getSAMRecord("TATT", "2M2D2M", q, true);

    final MachineErrorChooserInterface machineErrorChooserInterface = getMachineCycleCalibrator();
    final VariantAlignmentRecord r = new VariantAlignmentRecord(rec, 0, machineErrorChooserInterface, 0);
    assertTrue(Arrays.equals(new byte[] {0, 1, 2, 3}, r.getRecalibratedQuality()));
  }

  private MachineErrorChooserInterface getMachineCycleCalibrator() {
    return new MachineErrorChooserAdapter() {
        @Override
        public PhredScaler machineErrors(SAMReadGroupRecord rg, boolean readPaired) {
          return (qual, readPosition, arm) -> readPosition;
        }
      };
  }

  public void testRecalibrationReadArmLeft() {
    final byte[] q = {15, 21, 9, 10};
    final SAMRecord rec = getSAMRecord("TATT", "2M2D2M", q, true);

    final MachineErrorChooserInterface machineErrorChooserInterface = getArmBasedChooserInterface();
    final VariantAlignmentRecord r = new VariantAlignmentRecord(rec, 0, machineErrorChooserInterface, 0);
    assertTrue(Arrays.equals(new byte[] {30, 42, 18, 20}, r.getRecalibratedQuality()));
  }

  public void testRecalibrationReadArmRight() {
    final byte[] q = {15, 21, 9, 10};
    final SAMRecord rec = getSAMRecord("TATT", "2M2D2M", q, false);

    final MachineErrorChooserInterface machineErrorChooserInterface = getArmBasedChooserInterface();
    final VariantAlignmentRecord r = new VariantAlignmentRecord(rec, 0, machineErrorChooserInterface, 0);
    assertTrue(Arrays.equals(new byte[] {15, 21, 9, 10}, r.getRecalibratedQuality()));
  }

  private MachineErrorChooserInterface getArmBasedChooserInterface() {
    return new MachineErrorChooserAdapter() {
        @Override
        public PhredScaler machineErrors(SAMReadGroupRecord rg, boolean readPaired) {
          return (qual, readPosition, arm) -> arm == Arm.LEFT ? qual * 2 : qual;
        }
      };
  }


  private SAMRecord getSAMRecord(String readString, String cigar, byte[] q, boolean first) {
    assert readString.length() == q.length : readString.length() + "" + q.length ;
    final SAMFileHeader header = new SAMFileHeader();
    final SAMReadGroupRecord readGroup = new SAMReadGroupRecord(READ_GROUP_ID);
    readGroup.setPlatform(MachineType.ILLUMINA_PE.platform());
    final SAMReadGroupRecord cgReadGroup = new SAMReadGroupRecord(CG1_READ_GROUP_ID);
    cgReadGroup.setPlatform(MachineType.COMPLETE_GENOMICS.platform());
    final SAMReadGroupRecord cg2ReadGroup = new SAMReadGroupRecord(CG1_READ_GROUP_ID);
    cg2ReadGroup.setPlatform(MachineType.COMPLETE_GENOMICS.platform());
    header.addReadGroup(readGroup);
    header.addReadGroup(cgReadGroup);
    final SAMRecord rec = new SAMRecord(header);
    rec.setReadString(readString);
    rec.setReadPairedFlag(true);
    rec.setFirstOfPairFlag(first);
    rec.setSecondOfPairFlag(!first);
    rec.setAlignmentStart(42);
    rec.setCigarString(cigar);
    rec.setMappingQuality(43);
    rec.setAttribute("NH", 2);
    rec.setBaseQualities(q);
    rec.setInferredInsertSize(-232);
    rec.getReadGroup();
    rec.setAttribute(ReadGroupUtils.RG_ATTRIBUTE, READ_GROUP_ID);
    return rec;
  }

  private SAMRecord getCgSAMRecord(String readString, String cgOverlapBases, String cigar, String superCigar, byte[] q, String xq, boolean first, boolean reverse) {
    final SAMRecord samRecord = getSAMRecord(readString, cigar, q, first);
    samRecord.setAttribute(SamUtils.CG_SUPER_CIGAR, superCigar);
    samRecord.setAttribute(SamUtils.ATTRIBUTE_CG_OVERLAP_BASES, cgOverlapBases);
    samRecord.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, xq);
    samRecord.setAttribute(ReadGroupUtils.RG_ATTRIBUTE, CG2_READ_GROUP_ID);
    samRecord.setReadNegativeStrandFlag(reverse);
    return samRecord;
  }

  public void testCg1ForwardOverlap1() {
    final int overlap = 1;
    final int samReadLength = CgUtils.CG_RAW_READ_LENGTH - overlap;

    final byte[] q = new byte[samReadLength];
    Arrays.fill(q, (byte) 10);
    final String read = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
    final String overlapBases = "A";
    final String xq = FastaUtils.rawToAsciiString(new byte[]  {55});
    final SAMRecord rec = getCgSAMRecord(read, overlapBases, "24=7N10=", "5=1B20=7N10=", q, xq, true, false);

    final MachineErrorChooserInterface machineErrorChooserInterface = getMachineCycleCalibrator();
    final VariantAlignmentRecord r = new VariantAlignmentRecord(rec, 0, machineErrorChooserInterface, 0);
    final byte[] expected = new byte[samReadLength];
    for (int i = 0; i < expected.length; ++i) {
      expected[i] = (byte) (i + (i > 4 ? 1 : 0));
    }

    TestUtils.assertEquals(expected, r.getRecalibratedQuality());
    TestUtils.assertEquals(new byte[] {5}, r.getOverlapQuality());
  }

  public void testCg1ReverseOverlap1() {
    final int overlap = 1;
    final int samReadLength = CgUtils.CG_RAW_READ_LENGTH - overlap;

    final byte[] q = new byte[samReadLength];
    Arrays.fill(q, (byte) 10);
    final String read = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
    final String overlapBases = "A";
    final String xq = FastaUtils.rawToAsciiString(new byte[]  {55});
    final SAMRecord rec = getCgSAMRecord(read, overlapBases, "10=7N24=", "10=7N20=1B5=", q, xq, true, true);

    final MachineErrorChooserInterface machineErrorChooserInterface = getMachineCycleCalibrator();
    final VariantAlignmentRecord r = new VariantAlignmentRecord(rec, 0, machineErrorChooserInterface, 0);
    final byte[] expected = new byte[samReadLength];
    final int overlapPosition = samReadLength - 1 - CgUtils.CG_OVERLAP_POSITION;
    for (int i = 0; i < expected.length; ++i) {
      expected[i] = (byte) (CgUtils.CG_RAW_READ_LENGTH - i - (i > overlapPosition ? 1 : 0) - 1);
    }

    TestUtils.assertEquals(expected, r.getRecalibratedQuality());
    TestUtils.assertEquals(new byte[] {5}, r.getOverlapQuality());
  }


  public void testCg1ForwardOverlap3() {
    final int overlap = 3;
    final int samReadLength = CgUtils.CG_RAW_READ_LENGTH - overlap;
    final byte[] q = new byte[samReadLength];
    Arrays.fill(q, (byte) 10);
    final String read = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
    final String overlapBases = "AAA";
    final String xq = FastaUtils.rawToAsciiString(new byte[]  {55, 55, 55});
    final SAMRecord rec = getCgSAMRecord(read, overlapBases, "22=7N10=", "5=3B20=7N10=", q, xq, true, false);

    final MachineErrorChooserInterface machineErrorChooserInterface = getMachineCycleCalibrator();
    final VariantAlignmentRecord r = new VariantAlignmentRecord(rec, 0, machineErrorChooserInterface, 0);
    final byte[] expected = new byte[samReadLength];
    for (int i = 0; i < expected.length; ++i) {
      expected[i] = (byte) (i + (i > 4 ? 3 : 0));
    }

    TestUtils.assertEquals(expected, r.getRecalibratedQuality());
    TestUtils.assertEquals(new byte[] {5, 6, 7}, r.getOverlapQuality());
  }

  public void testCg1ReverseOverlap3() {
    final int overlap = 3;
    final int samReadLength = CgUtils.CG_RAW_READ_LENGTH - overlap;

    final byte[] q = new byte[samReadLength];
    Arrays.fill(q, (byte) 10);
    final String read = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
    final String overlapBases = "AAA";
    final String xq = FastaUtils.rawToAsciiString(new byte[]  {55, 55, 55});
    final SAMRecord rec = getCgSAMRecord(read, overlapBases, "10=7N22=", "10=7N20=2B5=", q, xq, true, true);

    final MachineErrorChooserInterface machineErrorChooserInterface = getMachineCycleCalibrator();
    final VariantAlignmentRecord r = new VariantAlignmentRecord(rec, 0, machineErrorChooserInterface, 0);
    final byte[] expected = new byte[samReadLength];
    final int overlapPosition = samReadLength - 1 - CgUtils.CG_OVERLAP_POSITION;
    for (int i = 0; i < expected.length; ++i) {
      expected[i] = (byte) (CgUtils.CG_RAW_READ_LENGTH - i - (i > overlapPosition ? 3 : 0) - 1);
    }

    TestUtils.assertEquals(expected, r.getRecalibratedQuality());
    TestUtils.assertEquals(new byte[] {7, 6, 5}, r.getOverlapQuality());
  }

  public void testCg1ForwardOverlap3SecondRead() {
    final int overlap = 3;
    final int samReadLength = CgUtils.CG_RAW_READ_LENGTH - overlap;

    final byte[] q = new byte[samReadLength];
    Arrays.fill(q, (byte) 10);
    final String read = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
    final String overlapBases = "AAA";
    final String xq = FastaUtils.rawToAsciiString(new byte[]  {55, 55, 55});
    final SAMRecord rec = getCgSAMRecord(read, overlapBases, "10=7N22=", "10=7N20=2B5=", q, xq, false, false);

    final MachineErrorChooserInterface machineErrorChooserInterface = getMachineCycleCalibrator();
    final VariantAlignmentRecord r = new VariantAlignmentRecord(rec, 0, machineErrorChooserInterface, 0);
    final byte[] expected = new byte[samReadLength];
    final int overlapPosition = samReadLength - 1 - CgUtils.CG_OVERLAP_POSITION;
    for (int i = 0; i < expected.length; ++i) {
      expected[i] = (byte) (i + (i > overlapPosition ? 3 : 0));
    }

    TestUtils.assertEquals(expected, r.getRecalibratedQuality());
    TestUtils.assertEquals(new byte[] {27, 28, 29}, r.getOverlapQuality());
  }

  public void testCg2ForwardOverlap2FirstRead() {
    final int overlap = 2;
    final int samReadLength = CgUtils.CG2_RAW_READ_LENGTH - overlap;

    final byte[] q = new byte[samReadLength];
    Arrays.fill(q, (byte) 10);
    final String read = "GGGGGGGGGGGGGGGGGGGGGGGGGGG";
    final String overlapBases = "AA";
    final String xq = FastaUtils.rawToAsciiString(new byte[]  {27, 28});
    final SAMRecord rec = getCgSAMRecord(read, overlapBases, "27=", "10=2B19=", q, xq, true, false);

    final MachineErrorChooserInterface machineErrorChooserInterface = getMachineCycleCalibrator();
    final VariantAlignmentRecord r = new VariantAlignmentRecord(rec, 0, machineErrorChooserInterface, 0);
    final byte[] expected = new byte[samReadLength];
    final int overlapPosition = CgUtils.CG2_OVERLAP_POSITION;
    for (int i = 0; i < expected.length; ++i) {
      expected[i] = (byte) (i + (i >= overlapPosition ? overlap : 0));
    }

    TestUtils.assertEquals(expected, r.getRecalibratedQuality());
    TestUtils.assertEquals(new byte[] {10, 11}, r.getOverlapQuality());
  }
  public void testCg2ForwardOverlap2SecondRead() {
    final int overlap = 2;
    final int samReadLength = CgUtils.CG2_RAW_READ_LENGTH - overlap;

    final byte[] q = new byte[samReadLength];
    Arrays.fill(q, (byte) 10);
    final String read = "GGGGGGGGGGGGGGGGGGGGGGGGGGG";
    final String overlapBases = "AA";
    final String xq = FastaUtils.rawToAsciiString(new byte[]  {27, 28});
    final SAMRecord rec = getCgSAMRecord(read, overlapBases, "27=", "10=2B19=", q, xq, false, false);

    final MachineErrorChooserInterface machineErrorChooserInterface = getMachineCycleCalibrator();
    final VariantAlignmentRecord r = new VariantAlignmentRecord(rec, 0, machineErrorChooserInterface, 0);
    final byte[] expected = new byte[samReadLength];
    final int overlapPosition = CgUtils.CG2_OVERLAP_POSITION;
    for (int i = 0; i < expected.length; ++i) {
      expected[i] = (byte) (i + (i >= overlapPosition ? overlap : 0));
    }

    TestUtils.assertEquals(expected, r.getRecalibratedQuality());
    TestUtils.assertEquals(new byte[] {10, 11}, r.getOverlapQuality());
  }

  public void testCg2ReverseOverlap2FirstRead() {
    final int overlap = 2;
    final int samReadLength = CgUtils.CG2_RAW_READ_LENGTH - overlap;

    final byte[] q = new byte[samReadLength];
    Arrays.fill(q, (byte) 10);
    final String read = "GGGGGGGGGGGGGGGGGGGGGGGGGGG";
    final String overlapBases = "AA";
    final String xq = FastaUtils.rawToAsciiString(new byte[]  {27, 28});
    final SAMRecord rec = getCgSAMRecord(read, overlapBases, "27=", "19=2B10=", q, xq, true, true);

    final MachineErrorChooserInterface machineErrorChooserInterface = getMachineCycleCalibrator();
    final VariantAlignmentRecord r = new VariantAlignmentRecord(rec, 0, machineErrorChooserInterface, 0);
    final byte[] expected = new byte[samReadLength];
    final int overlapPosition = samReadLength - CgUtils.CG2_OVERLAP_POSITION;
    for (int i = 0; i < expected.length; ++i) {
      expected[i] = (byte) (CgUtils.CG2_RAW_READ_LENGTH - 1 - (i + (i >= overlapPosition ? overlap : 0)));
    }

    TestUtils.assertEquals(expected, r.getRecalibratedQuality());
    TestUtils.assertEquals(new byte[] {11, 10}, r.getOverlapQuality());
  }
  public void testCg2ReverseOverlap2SecondRead() {
    final int overlap = 2;
    final int samReadLength = CgUtils.CG2_RAW_READ_LENGTH - overlap;

    final byte[] q = new byte[samReadLength];
    Arrays.fill(q, (byte) 10);
    final String read = "GGGGGGGGGGGGGGGGGGGGGGGGGGG";
    final String overlapBases = "AA";
    final String xq = FastaUtils.rawToAsciiString(new byte[]  {27, 28});
    final SAMRecord rec = getCgSAMRecord(read, overlapBases, "27=", "19=2B10=", q, xq, false, true);

    final MachineErrorChooserInterface machineErrorChooserInterface = getMachineCycleCalibrator();
    final VariantAlignmentRecord r = new VariantAlignmentRecord(rec, 0, machineErrorChooserInterface, 0);
    final byte[] expected = new byte[samReadLength];
    final int overlapPosition = samReadLength - CgUtils.CG2_OVERLAP_POSITION;
    for (int i = 0; i < expected.length; ++i) {
      expected[i] = (byte) (CgUtils.CG2_RAW_READ_LENGTH - 1 - (i + (i >= overlapPosition ? overlap : 0)));
    }

    TestUtils.assertEquals(expected, r.getRecalibratedQuality());
    TestUtils.assertEquals(new byte[] {11, 10}, r.getOverlapQuality());
  }
  public void testMinQuality() {
    final String read = "GGGGGGGGGGGGGGGGGGGGGGGGGGG";
    final byte[] q = new byte[read.length()];
    Arrays.fill(q, (byte) (FastaUtils.PHRED_LOWER_LIMIT_CHAR + 30));
    q[5] = 14;
    q[15] = 14;
    final SAMRecord rec = getSAMRecord(read, read.length() + "=", q, true);
    final VariantAlignmentRecord r = new VariantAlignmentRecord(rec, 0, new DefaultMachineErrorChooser(), 20);
    assertEquals('G', DnaUtils.getBase(r.getRead()[4]));
    assertEquals('N', DnaUtils.getBase(r.getRead()[5]));
    assertEquals('N', DnaUtils.getBase(r.getRead()[15]));
  }
}
