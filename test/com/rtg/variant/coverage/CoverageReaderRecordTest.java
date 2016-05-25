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

package com.rtg.variant.coverage;

import java.util.BitSet;

import com.rtg.sam.SamUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import junit.framework.TestCase;

/**
 */
public class CoverageReaderRecordTest extends TestCase {

  public void testStuff() {

    final SAMRecord blankSam = new SAMRecord(null);
    blankSam.setCigarString("1=");
    blankSam.setAttribute(SamUtils.ATTRIBUTE_IH, 2);
    final CoverageReaderRecord crr0 = new CoverageReaderRecord(blankSam, 0, false);
    assertEquals(-1, crr0.getFragmentLength());
    assertEquals(-1, crr0.getMateSequenceId());
    assertFalse(crr0.isMated());
    assertEquals(0.5, crr0.getCoverageMultiplier());
    crr0.setNextInChain(crr0);
    assertEquals(crr0, crr0.chain());
    assertEquals(0, crr0.compareTo(crr0));
    assertTrue(crr0.equals(crr0));

    final SAMFileHeader sfh = new SAMFileHeader();
    final SAMSequenceDictionary dict = new SAMSequenceDictionary();
    dict.addSequence(new SAMSequenceRecord("slkhr", 34));
    sfh.setSequenceDictionary(dict);
    final SAMRecord sam = new SAMRecord(sfh);
    sam.setCigarString("1=10H1S1P1I1N11=1D1M1=");
    sam.setReadPairedFlag(true);
    sam.setProperPairFlag(true);
    sam.setMateAlignmentStart(423);
    sam.setMateReferenceIndex(0);
    sam.setAlignmentStart(4);
    sam.setInferredInsertSize(460);

    final CoverageReaderRecord crr = new CoverageReaderRecord(sam, 23, false);
    assertEquals(23, crr.getGenome());
    assertEquals(460, crr.getFragmentLength());
    assertEquals(0, crr.getMateSequenceId());
    assertEquals(3, crr.getStart());
    assertEquals(16, crr.getLength());
    assertEquals(-1, crr.getSequenceId());
    assertEquals(1.0, crr.getCoverageMultiplier());
    BitSet coverage = crr.getCoverageBitSet();
    assertTrue(coverage.get(0));
    assertFalse(coverage.get(1));
    assertTrue(coverage.get(2));
    assertTrue(coverage.get(3));
    assertTrue(coverage.get(4));
    assertTrue(coverage.get(5));
    assertTrue(coverage.get(6));
    assertTrue(coverage.get(7));
    assertTrue(coverage.get(8));
    assertTrue(coverage.get(9));
    assertTrue(coverage.get(10));
    assertTrue(coverage.get(11));
    assertTrue(coverage.get(12));

    assertFalse(coverage.get(13));
    assertTrue(coverage.get(14));
    assertTrue(coverage.get(15));
    assertFalse(coverage.get(16));
    assertEquals(-1, crr0.compareTo(crr));
    assertFalse(crr0.equals(crr));

    coverage = new CoverageReaderRecord(sam, 23, true).getCoverageBitSet();
    assertTrue(coverage.get(0));
    assertFalse(coverage.get(1));
    assertTrue(coverage.get(2));
    assertTrue(coverage.get(3));
    assertTrue(coverage.get(4));
    assertTrue(coverage.get(5));
    assertTrue(coverage.get(6));
    assertTrue(coverage.get(7));
    assertTrue(coverage.get(8));
    assertTrue(coverage.get(9));
    assertTrue(coverage.get(10));
    assertTrue(coverage.get(11));
    assertTrue(coverage.get(12));

    assertTrue(coverage.get(13));
    assertTrue(coverage.get(14));
    assertTrue(coverage.get(15));
    assertFalse(coverage.get(16));

  }
}
