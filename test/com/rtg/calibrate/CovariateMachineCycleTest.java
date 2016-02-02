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
public class CovariateMachineCycleTest extends TestCase {

  public void testSam() throws BadSuperCigarException {
    final Covariate var = new CovariateMachineCycle(15);
    assertEquals(15, var.size());
    final SAMRecord sam = new SAMRecord(null);

    final Calibrator calwtf = new Calibrator(new Covariate[] {new CovariateReadGroup()}, null);
    final CalibratorCigarParser parser = new CalibratorCigarParser(calwtf) {
      @Override
      public int getReadPosition() {
        return 14;
      }
    };
    assertEquals(14, var.value(sam, parser));
    assertEquals("14", var.valueString(14));
  }

  public void testResize() throws BadSuperCigarException {
    final Covariate var = new CovariateMachineCycle(6);
    assertEquals(6, var.size());
    final SAMRecord sam = new SAMRecord(null);
    sam.setReadBases(new byte[7]);

    final Calibrator calwtf = new Calibrator(new Covariate[] {new CovariateReadGroup()}, null);
    final CalibratorCigarParser parser = new CalibratorCigarParser(calwtf) {
      private int mPos = 3;
      @Override
      public int getReadPosition() {
        mPos *= 2;
        return mPos;
      }
    };
    assertEquals(6, var.value(sam, parser));
    assertEquals(7, var.newSize());
    assertEquals(6, var.size());
    assertTrue(var.sizeChanged());
    sam.setReadBases(new byte[13]);
    assertEquals(12, var.value(sam, parser));
    assertEquals(13, var.newSize());
    assertEquals(6, var.size());
    assertTrue(var.sizeChanged());
    var.resized();
    assertEquals(13, var.size());
  }

  public void testNameType() {
    final Covariate var = new CovariateMachineCycle(15);
    assertEquals("machinecycle:15", var.name());
    assertEquals(CovariateEnum.MACHINECYCLE, var.getType());
  }
}
