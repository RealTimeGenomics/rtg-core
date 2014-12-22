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
package com.rtg.util.array.objectindex;


import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Test Create
 */
public class ObjectCreateTest extends TestCase {

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(ObjectCreateTest.class);
  }

  public static Test suite() {
    return new TestSuite(ObjectCreateTest.class);
  }

  /**
   * Constructor for CreateTest.
   */
  public ObjectCreateTest(final String arg0) {
    super(arg0);
  }

  public void testBad() {
    try {
      ObjectCreate.<Integer>createIndex(-1);
      fail("NegativeArraySizeException expected");
    } catch (final NegativeArraySizeException e) {
      //expected
      assertEquals("Negative length=-1", e.getMessage());
    }

    final ObjectIndex<Integer> index = ObjectCreate.createIndex(0);
    index.integrity();
    final ObjectIndex<Integer> i = ObjectCreate.createIndex(20);
    boolean test = i instanceof ObjectArray;
    assertTrue(test);
  }
}


