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
package com.rtg.simulation;

import junit.framework.TestCase;

/**
 */
public class ReadMappingRocTest extends TestCase {

  public void testROC() {
    final ReadMappingRoc roc = new ReadMappingRoc("test");
    assertEquals(0, roc.getTp(20), 0.001);
    assertEquals(0, roc.getMaxScore());

    roc.addTp(1);
    roc.addTp(1);
    roc.addTp(1);
    roc.addTp(1);
    roc.addFp(2);
    roc.addFp(2);
    roc.addFp(3);
    roc.addFp(3);

    assertEquals(3, roc.getMaxScore());

    assertEquals(0, roc.getTp(0), 0.001);
    assertEquals(0, roc.getFp(0), 0.001);

    assertEquals(4, roc.getTp(1), 0.001);
    assertEquals(0, roc.getFp(1), 0.001);

    assertEquals(0, roc.getTp(2), 0.001);
    assertEquals(2, roc.getFp(2), 0.001);

    assertEquals(0, roc.getTp(3), 0.001);
    assertEquals(2, roc.getFp(3), 0.001);

    roc.addFp(3, 0.25);
    assertEquals(2.25, roc.getFp(3), 0.001);

    //System.err.println(roc.getDistribution());
    //System.err.println(roc.getRoc());
  }

}
