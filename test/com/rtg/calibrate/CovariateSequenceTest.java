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
package com.rtg.calibrate;

import net.sf.samtools.SAMRecord;

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
