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


import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests the corresponding class.
 *
 */
public class RandomDnaTest extends TestCase {

  /**
   * Constructor (needed for JUnit)
   * @param name A string which names the object.
   */
  public RandomDnaTest(final String name) {
    super(name);
  }
  public static Test suite() {
    return new TestSuite(RandomDnaTest.class);
  }
  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  public void test() {
    final String s = RandomDna.random(10000);
    assertEquals(10000, s.length());
    int a = 0;
    int c = 0;
    int g = 0;
    int t = 0;
    for (int i = 0; i < 10000; i++) {
      switch (s.charAt(i)) {
      case 'A':
        a++;
        break;
      case 'C':
        c++;
        break;
      case 'G':
        g++;
        break;
      case 'T':
        t++;
        break;
      default:
        fail();
      }
    }
    assertTrue(a > 1000);
    assertTrue(c > 1000);
    assertTrue(g > 1000);
    assertTrue(t > 1000);
    assertTrue(a < 3500);
    assertTrue(c < 3500);
    assertTrue(g < 3500);
    assertTrue(t < 3500);
  }
}

