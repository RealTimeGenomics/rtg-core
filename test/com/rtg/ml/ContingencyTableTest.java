/*
 * Copyright (c) 2015. Real Time Genomics Limited.
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
public class ContingencyTableTest extends TestCase {

  public void test() {
    ContingencyTable eval = new ContingencyTable(100, 50, 25, 10);

    assertEquals(100.0, eval.truePositives());
    assertEquals(50.0, eval.falsePositives());
    assertEquals(25.0, eval.trueNegatives());
    assertEquals(10.0, eval.falseNegatives());
    assertEquals(125.0, eval.correct());
    assertEquals(60.0, eval.incorrect());

    assertEquals(185.0, eval.total());

    assertEquals(0.67, eval.accuracy(), 0.01);
    assertEquals(0.33, eval.errorRate(), 0.01);

    final double precision = ContingencyTable.precision(eval.truePositives(), eval.falsePositives());
    assertEquals(0.67, precision, 0.01);
    final double recall = ContingencyTable.recall(eval.truePositives(), eval.falseNegatives());
    assertEquals(0.91, recall, 0.01);
    assertEquals(0.77, ContingencyTable.fMeasure(precision, recall), 0.01);
  }
}
