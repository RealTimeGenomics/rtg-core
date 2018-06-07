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
package com.rtg.variant.match;

import java.io.IOException;
import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.rtg.reader.FastaUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.ReadGroupMachineErrorChooser;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.dna.DNARangeNAT;
import com.rtg.variant.util.VariantUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;


/**
 */
public class AlignmentMatchTest extends TestCase {

  private static final String TEST_READ_STRING = "ACGT"; // We use the "CGT" part
  private static final String TEST_READ_QUALITY_STRING = "!!%`"; // 0,0,4,63
  private static final byte[] TEST_RAW_QUALITY = FastaUtils.asciiToRawQuality(TEST_READ_QUALITY_STRING.toCharArray());

  public void test0() {
    final AlignmentMatch ins = new AlignmentMatch(null, "", "", 0, 0, 0, 7);
    ins.integrity();
    assertTrue(ins.isFixedLeft());
    assertTrue(ins.isFixedRight());

    assertEquals("", ins.toString());
    try {
      ins.read(0);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      // ok
    }
    try {
      ins.baseError(0);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      // ok
    }
    assertEquals(VariantUtils.phredToProb(7), ins.mapError(), 1E-14);
    assertEquals("", ins.readString());
    assertEquals(0, ins.length());
  }

  private SAMFileHeader getSAMFileHeader() {
    final SAMFileHeader hdr = new SAMFileHeader();
    final SAMReadGroupRecord illumina = new SAMReadGroupRecord("Illumina");
    illumina.setSample("sample1");
    illumina.setPlatform("illumina");
    hdr.addReadGroup(illumina);
    final SAMReadGroupRecord cg = new SAMReadGroupRecord("CG");
    cg.setSample("sample2");
    cg.setPlatform("complete");
    hdr.addReadGroup(cg);
    return hdr;
  }

  public void testNoQualities() {
    Diagnostic.setLogStream();
    final SAMFileHeader hdr = getSAMFileHeader();
    final SAMRecord rec = new SAMRecord(hdr);
    final VariantAlignmentRecord var = new VariantAlignmentRecord(rec);


    // no qualities
    final AlignmentMatch ins = new AlignmentMatch(var, null, TEST_READ_STRING, null, 20, 1, 3, 7, true, true);
    ins.integrity();
    assertEquals("CGT", ins.toString());
    assertEquals(3, ins.length());
    assertEquals(DNARangeNAT.C, ins.read(0));
    assertEquals(DNARangeNAT.T, ins.read(2));
    assertEquals(VariantUtils.phredToProb(20), ins.baseError(0), 1E-14);
    assertEquals(VariantUtils.phredToProb(20), ins.baseError(1), 1E-14);
    assertEquals(VariantUtils.phredToProb(20), ins.baseError(2), 1E-14);
    assertEquals(VariantUtils.phredToProb(7), ins.mapError(), 1E-14);
  }

  public void testQualities() {
    Diagnostic.setLogStream();
    final SAMFileHeader hdr = getSAMFileHeader();
    final SAMRecord rec = new SAMRecord(hdr);
    final VariantAlignmentRecord var = new VariantAlignmentRecord(rec);

    final AlignmentMatch ins = new AlignmentMatch(var, null, TEST_READ_STRING, TEST_RAW_QUALITY, 20, 1, 3, 7, true, true);
    assertEquals(VariantUtils.phredToProb(0), ins.baseError(0), 1E-14);
    assertEquals(VariantUtils.phredToProb(4), ins.baseError(1), 1E-14);
    assertEquals(VariantUtils.phredToProb(63), ins.baseError(2), 1E-14);
    assertEquals(VariantUtils.phredToProb(7), ins.mapError(), 1E-14);
  }

  public void testNoReadGroupInSam() throws IOException {
    Diagnostic.setLogStream();
    final SAMFileHeader hdr = getSAMFileHeader();
    final SAMRecord rec = new SAMRecord(hdr);
    final VariantAlignmentRecord var = new VariantAlignmentRecord(rec);
    final ReadGroupMachineErrorChooser chooser = new ReadGroupMachineErrorChooser(hdr);

    // no read group in SAM record, so no quality correction curve
    try {
      new AlignmentMatch(var, chooser, TEST_READ_STRING, TEST_RAW_QUALITY, 20, 1, 3, 7, true, true);
      fail();
    } catch (final NoTalkbackSlimException ntse) {
      //expected
    }
  }

  public void testIlluminaCorrectionCurve() throws IOException {
    Diagnostic.setLogStream();
    final SAMFileHeader hdr = getSAMFileHeader();
    final SAMRecord rec = new SAMRecord(hdr);
    final ReadGroupMachineErrorChooser chooser = new ReadGroupMachineErrorChooser(hdr);

    // Illumina correction curve
    rec.setAttribute("RG", "Illumina");
    final AlignmentMatch ins = new AlignmentMatch(new VariantAlignmentRecord(rec), chooser, TEST_READ_STRING, TEST_RAW_QUALITY, 20, 1, 3, 7, true, true);
    ins.integrity();
    assertEquals(VariantUtils.phredToProb(0), ins.baseError(0), 1E-14);
    assertEquals(VariantUtils.phredToProb(4), ins.baseError(1), 1E-14);
    assertEquals(VariantUtils.phredToProb(63), ins.baseError(2), 1E-14);
    assertEquals(VariantUtils.phredToProb(7), ins.mapError(), 1E-14);

  }

  public void testCgCorrectionCurve() throws IOException {
    Diagnostic.setLogStream();
    final SAMFileHeader hdr = getSAMFileHeader();
    final SAMRecord rec = new SAMRecord(hdr);
    final ReadGroupMachineErrorChooser chooser = new ReadGroupMachineErrorChooser(hdr);

    // Complete Genomics correction curve
    rec.setAttribute("RG", "CG");
    rec.setBaseQualities(TEST_RAW_QUALITY);
    final VariantAlignmentRecord var = new VariantAlignmentRecord(rec, 0, chooser, 0);
    final AlignmentMatch ins = new AlignmentMatch(var, chooser, TEST_READ_STRING, var.getRecalibratedQuality(), 20, 1, 3, 7, true, true);
    ins.integrity();
    assertEquals(VariantUtils.phredToProb(0), ins.baseError(0), 1E-14);
    assertEquals(VariantUtils.phredToProb(5), ins.baseError(1), 1E-14);
    assertEquals(VariantUtils.phredToProb(10), ins.baseError(2), 1E-14);
    assertEquals(VariantUtils.phredToProb(7), ins.mapError(), 1E-14);
  }

  public void testReadRangeLower() {
    final AlignmentMatch ins = getTestAlignmentMatch();
    try {
      ins.read(-1);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      // ok
    }
  }
  public void testReadRangeUpper() {
    final AlignmentMatch ins = getTestAlignmentMatch();
    try {
      ins.read(3);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      // ok
    }
  }

  public void testReadValues() {
    final AlignmentMatch ins = getTestAlignmentMatch();
    assertEquals(Arrays.asList(3, 0, 1),
      IntStream.range(0, 3).map(ins::read).boxed().collect(Collectors.toList())
    );
  }

  public void testErrorRangeLower() {
    final AlignmentMatch ins = getTestAlignmentMatch();
    try {
      ins.baseError(-1);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      // ok
    }
  }

  public void testErrorRangeUpper() {
    final AlignmentMatch ins = getTestAlignmentMatch();
    try {
      ins.baseError(3);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      // ok
    }
  }

  public void testBaseError() {
    final AlignmentMatch ins = getTestAlignmentMatch();
    assertEquals(VariantUtils.phredToProb('D' - FastaUtils.PHRED_LOWER_LIMIT_CHAR), ins.baseError(0), 1E-14);
    assertEquals(VariantUtils.phredToProb('A' - FastaUtils.PHRED_LOWER_LIMIT_CHAR), ins.baseError(1), 1E-14);
    assertEquals(VariantUtils.phredToProb('B' - FastaUtils.PHRED_LOWER_LIMIT_CHAR), ins.baseError(2), 1E-14);
  }

  public void testFixedEnds() {
    final AlignmentMatch ins = getTestAlignmentMatch();
    assertFalse(ins.isFixedLeft());
    assertTrue(ins.isFixedRight());

  }
  public void testToString() {
    assertEquals("~TAC", getTestAlignmentMatch().toString());
  }

  public void testMapError() {
    final AlignmentMatch ins = getTestAlignmentMatch();
    assertEquals(VariantUtils.phredToProb(7), ins.mapError(), 1E-14);
  }

  public void testIntegrity() {
    final AlignmentMatch ins = getTestAlignmentMatch();
    ins.integrity();
  }

  public void testReadString() {
    final AlignmentMatch ins = getTestAlignmentMatch();
    assertEquals("TAC", ins.readString());
  }
  public void testLength() {
    final AlignmentMatch ins = getTestAlignmentMatch();
    assertEquals(3, ins.length());
  }


  private AlignmentMatch getTestAlignmentMatch() {
    return new AlignmentMatch(null, null, "ACGTACGT", FastaUtils.asciiToRawQuality("ABCDABCD"), 0, 3, 3, 7, false, true);
  }
}

