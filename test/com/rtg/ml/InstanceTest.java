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
package com.rtg.ml;

import junit.framework.TestCase;

/**
 */
public class InstanceTest extends TestCase {

  public void test() {
    final double[] data = {};
    final Instance i = new Instance(data, true);
    assertTrue(data == i.instance());
    assertTrue(i.isPositive());
    final Instance j = i.reweight(0.5);
    assertEquals(0.5, j.weight(), 1e-10);
  }

}
