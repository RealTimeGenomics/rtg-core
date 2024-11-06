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

import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class CovariateSequenceTest extends TestCase {

  private SAMRecord recSam(String sequence) {
    final SAMRecord sam = new SAMRecord(null);
    sam.setReferenceName(sequence);
    sam.setCigarString("35=");
    sam.setReadString("acgtacgtggacgtacgtggacgtacgtggttttt");
    sam.setBaseQualityString("012345678901234567890123456789ABCDE");
    sam.setAlignmentStart(1);
    sam.setMappingQuality(1);
    return sam;
  }

  public void testNewSize() {
    final CovariateSequence cs = new CovariateSequence();
    final Calibrator cal = new Calibrator(new Covariate[] {cs, new CovariateReadGroup()}, null);
    final CalibratorCigarParser parser = new CalibratorCigarParser(cal);
    cs.value(recSam("sequence1"), parser);
    assertEquals(1, cs.newSize());
    cs.value(recSam("sequence1"), parser);
    assertEquals(1, cs.newSize());
    cs.value(recSam("sequence2"), parser);
    assertEquals(2, cs.newSize());
    cs.value(recSam("sequence1"), parser);
    assertEquals(2, cs.newSize());
  }

  public void testValue() {
    final CovariateSequence cs = new CovariateSequence();
    final Calibrator cal = new Calibrator(new Covariate[] {cs, new CovariateReadGroup()}, null);
    final CalibratorCigarParser parser = new CalibratorCigarParser(cal);

    assertEquals(0, cs.value(recSam("sequence1"), parser));
    assertEquals(1, cs.value(recSam("sequence2"), parser));
    assertEquals(2, cs.value(recSam("sequence3"), parser));
    assertEquals(0, cs.value(recSam("sequence1"), parser));
    assertEquals(2, cs.value(recSam("sequence3"), parser));

  }

  public void testValueString() {
    final CovariateSequence cs = new CovariateSequence();
    final Calibrator cal = new Calibrator(new Covariate[] {cs, new CovariateReadGroup()}, null);
    final CalibratorCigarParser parser = new CalibratorCigarParser(cal);

    cs.value(recSam("sequence1"), parser);
    cs.value(recSam("sequence2"), parser);
    cs.value(recSam("sequence3"), parser);
    assertEquals("sequence1", cs.valueString(0));
    assertEquals("sequence2", cs.valueString(1));
    assertEquals("sequence3", cs.valueString(2));

  }

  public void testParse() {
    final CovariateSequence cs = new CovariateSequence();
    assertEquals(0, cs.parse("sequence1"));
    assertEquals(1, cs.parse("sequence2"));
    assertEquals(2, cs.parse("sequence3"));
    assertEquals(0, cs.parse("sequence1"));
    assertEquals(2, cs.parse("sequence3"));

  }

  public void testGetType() {
    final CovariateSequence cs = new CovariateSequence();
    assertEquals(CovariateEnum.SEQUENCE, cs.getType());
  }
}
