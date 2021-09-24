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
public class SimpleEvaluationTest extends TestCase {

  public void test() {
    final SimpleEvaluation eval = new SimpleEvaluation();

    final BuildClassifier b = new ZeroRBuilder();
    final Dataset data = TrainTestSplitTest.makeSimpleDataset(100, 200);
    b.build(data);

    final PredictClassifier p = b.getClassifier();
    eval.evaluate(p, data);
    assertEquals(200.0, eval.correct());
    assertEquals(200.0, eval.trueNegatives());
    assertEquals(100.0, eval.incorrect());
    assertEquals(100.0, eval.falseNegatives());

    assertEquals(0.0, eval.truePositives());
    assertEquals(0.0, eval.falsePositives());

    assertEquals(300.0, eval.total());

    assertEquals(0.67, eval.accuracy(), 0.01);
    assertEquals(0.33, eval.errorRate(), 0.01);
  }
}
