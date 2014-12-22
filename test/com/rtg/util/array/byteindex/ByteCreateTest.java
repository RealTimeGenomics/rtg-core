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
package com.rtg.util.array.byteindex;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Test Create
 */
public class ByteCreateTest extends TestCase {

  private static final long FREE_LIMIT = 2L * Integer.MAX_VALUE + 200000000L; //allow 200m of freeboard

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(ByteCreateTest.class);
  }

  public static Test suite() {
    return new TestSuite(ByteCreateTest.class);
  }

  /**
   * Constructor for CreateTest.
   */
  public ByteCreateTest(final String arg0) {
    super(arg0);
  }

  public void testBad() {
    try {
      ByteCreate.createIndex(-1);
      fail("NegativeArraySizeException expected");
    } catch (final NegativeArraySizeException e) {
      assertEquals("Negative length=-1", e.getMessage()); //expected
    }

    final ByteIndex index = ByteCreate.createIndex(0);
    index.integrity();
  }

  public void test() {
    final ByteIndex a = ByteCreate.createIndex(10);
    assertEquals(10, a.length());
    assertTrue(a instanceof ByteArray);

    //Only run if there is enough memory not to nuke everything
    System.gc();
    final long mem = Runtime.getRuntime().maxMemory();
    //System.err.println("FREE_LIMIT=" + FREE_LIMIT + " mem=" + mem);
    if (mem > FREE_LIMIT) {
      final ByteIndex b = ByteCreate.createIndex(ByteIndex.MAX_LENGTH);
      assertEquals(ByteIndex.MAX_LENGTH, b.length());
      assertTrue(b instanceof ByteArray);

      final ByteIndex c = ByteCreate.createIndex(ByteIndex.MAX_LENGTH + 1L);
      assertEquals(ByteIndex.MAX_LENGTH + 1L, c.length());
      assertTrue(c instanceof ByteChunks);
    }
    System.gc();
  }
}


