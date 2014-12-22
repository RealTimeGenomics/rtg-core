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
package com.rtg.util.array.shortindex;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Test Create
 */
public class ShortCreateTest extends TestCase {

  private static final long FREE_LIMIT = 2L * Integer.MAX_VALUE + 200000000L; //allow 200m of freeboard

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(ShortCreateTest.class);
  }

  public static Test suite() {
    return new TestSuite(ShortCreateTest.class);
  }

  /**
   * Constructor for CreateTest.
   */
  public ShortCreateTest(final String arg0) {
    super(arg0);
  }

  public void testBad() {
    try {
      ShortCreate.createIndex(-1);
      fail("NegativeArraySizeException expected");
    } catch (final NegativeArraySizeException e) {
      assertEquals("Negative length=-1", e.getMessage()); //expected
    }

    final ShortIndex index = ShortCreate.createIndex(0);
    index.integrity();
  }

  public void test() {
    final ShortIndex a = ShortCreate.createIndex(10);
    assertEquals(10, a.length());
    assertTrue(a instanceof ShortArray);

    //Only run if there is enough memory not to nuke everything
    System.gc();
    final long mem = Runtime.getRuntime().maxMemory();
    //System.err.println("FREE_LIMIT=" + FREE_LIMIT + " mem=" + mem);
    if (mem > FREE_LIMIT) {
      final ShortIndex b = ShortCreate.createIndex(ShortIndex.MAX_LENGTH);
      assertEquals(ShortIndex.MAX_LENGTH, b.length());
      assertTrue(b instanceof ShortArray);

      final ShortIndex c = ShortCreate.createIndex(ShortIndex.MAX_LENGTH + 1L);
      assertEquals(ShortIndex.MAX_LENGTH + 1L, c.length());
      assertTrue(c instanceof ShortChunks);
    }
    System.gc();
  }
}


