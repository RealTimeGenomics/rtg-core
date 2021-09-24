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

package com.rtg.variant.sv;

import junit.framework.TestCase;

/**
 */
public class SignalConstantLnTest extends TestCase {

  public void test() {
    final SamCounts sa = new SamArray(5);
    sa.increment(1);
    sa.increment(2);
    sa.increment(2);
    sa.increment(3);
    sa.increment(3);
    sa.increment(3);
    sa.increment(4);
    sa.increment(4);
    sa.increment(4);
    sa.increment(4);
    final SignalConstantLn sig = new SignalConstantLn(sa, new DistributionConstant(-2, 2, 2.0), "blah");
    final SignalDistributionLn sig2 = new SignalDistributionLn(sa, new DistributionConstant(-2, 2, 2.0), "blah");
    sig.globalIntegrity();

    assertEquals(6.3069, sig2.value(0), 0.0001);
    assertEquals(6.3069, sig.value(0), 0.0001);
    assertEquals(4.3069, sig2.value(1), 0.0001);
    assertEquals(4.3069, sig.value(1), 0.0001);
    assertEquals(2.5232, sig.value(2), 0.0001);
    assertEquals(1.2958, sig.value(3), 0.0001);
    assertEquals(2.9890, sig.value(4), 0.0001);
    assertEquals("blah", sig.columnLabel());
    assertEquals("blah", sig2.columnLabel());
  }
}
