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
package com.rtg.util.test;

import junit.framework.TestCase;

/**
 * Tests for the class with the according name
 */
public class NotRandomRandomTest extends TestCase {

  public NotRandomRandomTest(final String name) {
    super(name);
  }


  public void test() {
    NotRandomRandom rand = new NotRandomRandom();
    double dd = 0.0;
    for (int i = 0; i < 11; i++) {
      double r = rand.nextDouble();
      assertTrue(r > dd - 0.00000000001 && r < dd + 0.00000000001);
      dd += 0.1;
    }
    assertEquals(0.0, rand.nextDouble());
    assertEquals(0, rand.nextInt(3));
    assertEquals(1, rand.nextInt(3));
    assertEquals(2, rand.nextInt(3));
    assertEquals(0, rand.nextInt(3));
  }

}
