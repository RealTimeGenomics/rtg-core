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
package com.rtg.util;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 */
public class PairTest extends TestCase {

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(PairTest.class);
  }

  public static Test suite() {
    return new TestSuite(PairTest.class);
  }

  /*
   * Test method for 'com.reeltwo.art.generator.Pair'
   */
  public final void test() {
    final Pair<String, Integer> p = new Pair<>("abc", 0);
    assertEquals("abc:0", p.toString());
    assertEquals("abc", p.getA());
    assertEquals(Integer.valueOf(0), p.getB());
    assertEquals(p, p);
    assertFalse(p.equals(null));
    assertFalse(p.equals("abc"));

    final Pair<String, Integer> q = new Pair<>("abc", 0);
    assertEquals(p, q);
    assertEquals(p.hashCode(), q.hashCode());
    assertEquals(1388221096, p.hashCode()); // a regression test - hard to work out what it will be
  }

  /*
   * Test method for 'com.reeltwo.art.generator.Pair'
   */
  public final void testEquals() {
    final Object[][] groups = {
      {new Pair<>("", 0), new Pair<>("", 0) },
      {new Pair<>("", 1), new Pair<>("", 1) },
      {new Pair<>("a", 0), new Pair<>("a", 0) },
      {new Pair<>("a", 1), new Pair<>("a", 1) },
    };
    TestUtils.equalsHashTest(groups);
  }

  public final void testBad() {

    try {
      new Pair<String, String>(null, "");
      fail("NullPointerException expected");
    } catch (final NullPointerException e) {
      //expeted
    }
    try {
      new Pair<String, String>("", null);
      fail("NullPointerException expected");
    } catch (final NullPointerException e) {
      //expeted
    }
  }
}

