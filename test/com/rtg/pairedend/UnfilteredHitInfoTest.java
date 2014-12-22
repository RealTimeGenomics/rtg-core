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
package com.rtg.pairedend;


import junit.framework.TestCase;


/**
 */
public class UnfilteredHitInfoTest extends TestCase {

  public void testSetters() {
    final UnfilteredHitInfo uhi = new UnfilteredHitInfo();
    assertFalse(uhi.getMatedOk());
    assertFalse(uhi.getUnmatedOk());
    uhi.setMatedOk(true);
    uhi.setUnmatedOk(true);
    assertTrue(uhi.getMatedOk());
    assertTrue(uhi.getUnmatedOk());

    uhi.setValues(true, false, 0, 0);
    assertFalse(uhi.getMatedOk());
    assertFalse(uhi.getUnmatedOk());
  }
}
