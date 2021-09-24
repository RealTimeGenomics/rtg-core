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

import com.rtg.sam.BadSuperCigarException;

import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class CovariateBaseQualityTest extends TestCase {


  public void test() throws BadSuperCigarException {
    final Covariate var = new CovariateBaseQuality();
    assertEquals(64, var.size());
    assertEquals(64, var.newSize());
    final SAMRecord sam = new SAMRecord(null);
    sam.setBaseQualityString("ABCDEFGHIJ");

    final Calibrator calwtf = new Calibrator(new Covariate[] {new CovariateReadGroup()}, null);
    final CalibratorCigarParser parser = new CalibratorCigarParser(calwtf) {
      @Override
      public int getCurrentQuality() {
        return 'J' - '!';
      }
    };
    assertEquals('J' - '!', var.value(sam, parser));
    assertEquals("14", var.valueString(14));
    assertEquals(14, var.parse("14"));
  }

  public void testNoQuality() throws BadSuperCigarException {
    final Covariate var = new CovariateBaseQuality();
    final SAMRecord sam = new SAMRecord(null);
    sam.setBaseQualityString("*");

    final Calibrator calwtf = new Calibrator(new Covariate[] {new CovariateReadGroup()}, null);
    final CalibratorCigarParser parser = new CalibratorCigarParser(calwtf) {
      @Override
      public int getReadPosition() {
        return 0;
      }
    };
    assertEquals('5' - '!', var.value(sam, parser));
  }

  public void testNameType() {
    final Covariate var = new CovariateBaseQuality();
    assertEquals("basequality", var.name());
    assertEquals(CovariateEnum.BASEQUALITY, var.getType());
  }
}
