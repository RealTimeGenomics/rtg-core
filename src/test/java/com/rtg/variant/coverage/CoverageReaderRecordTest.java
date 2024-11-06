/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

  public void testEmptyCigar() {
    final SAMRecord sam = new SAMRecord(null);
    sam.setCigarString("*");
    sam.setAttribute(SamUtils.ATTRIBUTE_IH, 1);
    final CoverageReaderRecord record = new CoverageReaderRecord(sam, 0, false);
    assertEquals(-1, record.getFragmentLength());
    assertEquals(-1, record.getMateSequenceId());
    assertFalse(record.isMated());
    assertEquals(1.0, record.getCoverageMultiplier());
  }
}
