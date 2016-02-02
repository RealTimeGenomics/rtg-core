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
public class CovariateArmTest extends TestCase {

  public void testSamPaired() throws BadSuperCigarException {
    final Covariate var = new CovariateArm();
    assertEquals("arm", var.name());
    assertEquals(2, var.size());
    final SAMRecord sam = new SAMRecord(null);
    sam.setReadPairedFlag(true);
    sam.setFirstOfPairFlag(true);
    assertEquals(0, var.value(sam, null));
    sam.setFirstOfPairFlag(false);
    sam.setSecondOfPairFlag(true);
    assertEquals(1, var.value(sam, null));
  }
}
