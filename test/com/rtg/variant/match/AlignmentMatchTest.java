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

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.ReadGroupMachineErrorChooser;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.dna.DNARangeNAT;
import com.rtg.variant.util.VariantUtils;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;

import junit.framework.TestCase;


/**
 */
public class AlignmentMatchTest extends TestCase {

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

  public void testReadGroups() throws IOException {
    Diagnostic.setLogStream();
    final String read = "ACGT";   // we use the "CGT" part
    final String qualities = "!!%`";  // that is, 0,0,4,63.
    final SAMFileHeader hdr = new SAMFileHeader();
    final SAMReadGroupRecord illumina = new SAMReadGroupRecord("Illumina");
    illumina.setSample("sample1");
    illumina.setPlatform("illumina");
    hdr.addReadGroup(illumina);
    final SAMReadGroupRecord cg = new SAMReadGroupRecord("CG");
    cg.setSample("sample2");
    cg.setPlatform("complete");
    hdr.addReadGroup(cg);
    final SAMRecord rec = new SAMRecord(hdr);
    final VariantAlignmentRecord var = new VariantAlignmentRecord(rec);
    final ReadGroupMachineErrorChooser chooser = new ReadGroupMachineErrorChooser(hdr);

    // no qualities
    AlignmentMatch ins = new AlignmentMatch(var, null, read, null, 20, 1, 3, 7);
    ins.integrity();
    assertEquals("CGT", ins.toString());
    assertEquals(3, ins.length());
    assertEquals(DNARangeNAT.C, ins.read(0));
    assertEquals(DNARangeNAT.T, ins.read(2));
    assertEquals(VariantUtils.phredToProb(20), ins.baseError(0), 1E-14);
    assertEquals(VariantUtils.phredToProb(20), ins.baseError(1), 1E-14);
    assertEquals(VariantUtils.phredToProb(20), ins.baseError(2), 1E-14);
    assertEquals(VariantUtils.phredToProb(7), ins.mapError(), 1E-14);

    ins = new AlignmentMatch(var, null, read, qualities, 20, 1, 3, 7);
    assertEquals(VariantUtils.phredToProb(0), ins.baseError(0), 1E-14);
    assertEquals(VariantUtils.phredToProb(4), ins.baseError(1), 1E-14);
    assertEquals(VariantUtils.phredToProb(63), ins.baseError(2), 1E-14);
    assertEquals(VariantUtils.phredToProb(7), ins.mapError(), 1E-14);

    // no read group in SAM record, so no quality correction curve
    try {
      new AlignmentMatch(var, chooser, read, qualities, 20, 1, 3, 7);
      fail();
    } catch (final NoTalkbackSlimException ntse) {
      //expected
    }

    // Illumina correction curve
    rec.setAttribute("RG", "Illumina");
    ins = new AlignmentMatch(new VariantAlignmentRecord(rec), chooser, read, qualities, 20, 1, 3, 7);
    ins.integrity();
    assertEquals(VariantUtils.phredToProb(0), ins.baseError(0), 1E-14);
    assertEquals(VariantUtils.phredToProb(4), ins.baseError(1), 1E-14);
    assertEquals(VariantUtils.phredToProb(63), ins.baseError(2), 1E-14);
    assertEquals(VariantUtils.phredToProb(7), ins.mapError(), 1E-14);

    // Complete Genomics correction curve
    rec.setAttribute("RG", "CG");
    ins = new AlignmentMatch(new VariantAlignmentRecord(rec), chooser, read, qualities, 20, 1, 3, 7);
    ins.integrity();
    assertEquals(VariantUtils.phredToProb(0), ins.baseError(0), 1E-14);
    assertEquals(VariantUtils.phredToProb(5), ins.baseError(1), 1E-14);
    assertEquals(VariantUtils.phredToProb(10), ins.baseError(2), 1E-14);
    assertEquals(VariantUtils.phredToProb(7), ins.mapError(), 1E-14);
  }

  public void test() {
    final AlignmentMatch ins = new AlignmentMatch(null, "ACGTACGT", "ABCDABCD", 0, 3, 3, 7, false, true);
    ins.integrity();
    assertFalse(ins.isFixedLeft());
    assertTrue(ins.isFixedRight());

    try {
      ins.read(-1);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      // ok
    }
    try {
      ins.read(3);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      // ok
    }
    try {
      ins.baseError(-1);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      // ok
    }
    try {
      ins.baseError(4);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      // ok
    }

    assertEquals("~TAC", ins.toString());
    assertEquals(3, ins.read(0));
    assertEquals(0, ins.read(1));
    assertEquals(1, ins.read(2));
    assertEquals(VariantUtils.phredToProb('D' - '!'), ins.baseError(0), 1E-14);
    assertEquals(VariantUtils.phredToProb('A' - '!'), ins.baseError(1), 1E-14);
    assertEquals(VariantUtils.phredToProb('B' - '!'), ins.baseError(2), 1E-14);
    assertEquals(VariantUtils.phredToProb(7), ins.mapError(), 1E-14);
    assertEquals("TAC", ins.readString());
    assertEquals("~TAC", ins.toString());
    assertEquals(3, ins.length());
  }
}

