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
package com.rtg.sam;


import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests corresponding class
 */
public class PileUpTest extends TestCase {

  public void testPileUp() {


    PileUp pu = new PileUp(10);

    pu.add('t', 3);
    pu.add('a', 5);
    pu.add('n', 1);
    pu.add('n', 9);
    pu.add('a', 0);
    pu.add('g', 7);
    pu.add('g', 4);
    pu.add('c', 7);
    pu.add('c', 8);

    assertEquals(9, pu.total());
    assertEquals(0.9, pu.coverage());
    assertEquals(8, pu.consensus());
  }

  /**
   */
  public PileUpTest(final String name) {
    super(name);
  }
  public static Test suite() {
    return new TestSuite(PileUpTest.class);
  }
  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }


}
