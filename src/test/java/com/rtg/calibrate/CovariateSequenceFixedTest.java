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
package com.rtg.calibrate;

import java.util.ArrayList;

import com.rtg.sam.SamUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import junit.framework.TestCase;

/**
 */
public class CovariateSequenceFixedTest extends TestCase {

  private SAMRecord recSam(String sequence, SAMFileHeader header) {
    final SAMRecord sam = new SAMRecord(header);
    sam.setReferenceName(sequence);
    sam.setCigarString("35=");
    sam.setReadString("acgtacgtggacgtacgtggacgtacgtggttttt");
    sam.setBaseQualityString("012345678901234567890123456789ABCDE");
    sam.setAlignmentStart(1);
    sam.setMappingQuality(1);
    return sam;
  }

  private SAMFileHeader makeHeader() {
    final SAMFileHeader header = new SAMFileHeader();
    final SAMSequenceDictionary dict = header.getSequenceDictionary();
    dict.addSequence(new SAMSequenceRecord("sequence1", 100));
    dict.addSequence(new SAMSequenceRecord("sequence2", 100));
    dict.addSequence(new SAMSequenceRecord("sequence3", 100));
    return header;
  }

  public void testValue() {
    final SAMFileHeader header = makeHeader();
    final CovariateSequenceFixed cs = new CovariateSequenceFixed(SamUtils.getSequenceNames(header));
    final Calibrator cal = new Calibrator(new Covariate[] {cs, new CovariateReadGroup()}, null);
    final CalibratorCigarParser parser = new CalibratorCigarParser(cal);
    assertEquals(3, cs.newSize());
    assertEquals(3, cs.size());

    assertEquals(0, cs.value(recSam("sequence1", header), parser));
    assertEquals(1, cs.value(recSam("sequence2", header), parser));
    assertEquals(2, cs.value(recSam("sequence3", header), parser));
    assertEquals(0, cs.value(recSam("sequence1", header), parser));
    assertEquals(2, cs.value(recSam("sequence3", header), parser));

  }

  public void testValueString() {
    final SAMFileHeader header = makeHeader();
    final CovariateSequenceFixed cs = new CovariateSequenceFixed(SamUtils.getSequenceNames(header));
    final Calibrator cal = new Calibrator(new Covariate[] {cs, new CovariateReadGroup()}, null);
    final CalibratorCigarParser parser = new CalibratorCigarParser(cal);

    cs.value(recSam("sequence1", header), parser);
    cs.value(recSam("sequence2", header), parser);
    cs.value(recSam("sequence3", header), parser);
    assertEquals("sequence1", cs.valueString(0));
    assertEquals("sequence2", cs.valueString(1));
    assertEquals("sequence3", cs.valueString(2));

  }

  public void testParse() {
    final SAMFileHeader header = makeHeader();
    final CovariateSequenceFixed cs = new CovariateSequenceFixed(SamUtils.getSequenceNames(header));
    assertEquals(0, cs.parse("sequence1"));
    try {
      cs.parse("no-such-sequence1");
      fail();
    } catch (final UnsupportedOperationException e) {
      // ok
    }
  }

  public void testGetType() {
    final CovariateSequenceFixed cs = new CovariateSequenceFixed(new ArrayList<>());
    assertEquals(CovariateEnum.SEQUENCE, cs.getType());
  }
}
