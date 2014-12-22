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

import com.rtg.util.PortableRandom;

/**
 */
public class RandomTreeBuilderTest extends AbstractBuildClassifierTest {

  @Override
  BuildClassifier makeClassifier() {
    return new RandomTreeBuilder();
  }

  public void testSplitPoint() {
    assertEquals(10, (int) RandomTreeBuilder.getSplitPoint(MlDataType.INTEGER, 9, 11));
    assertEquals(10, (int) RandomTreeBuilder.getSplitPoint(MlDataType.INTEGER, 10, 11));

    assertEquals(10.0, RandomTreeBuilder.getSplitPoint(MlDataType.DOUBLE, 9.0, 11.0));
    assertEquals(10.5, RandomTreeBuilder.getSplitPoint(MlDataType.DOUBLE, 10.0, 11.0));

  }

  public void testTree() {
    BuildClassifier b = makeClassifier();
    Dataset data = TrainTestSplitTest.makeSimpleDataset(200, 100);

    TrainTestSplit split = TrainTestSplit.sampleWithReplacement(data, 150, new PortableRandom(42));
    b.build(split.mTrain);

    PredictClassifier p = b.getClassifier();

    String thetree = p.toString(new StringBuilder(), "", data).toString();
    assertNotNull(thetree);
    //System.err.println(thetree);
  }

  public void testCircleTree() {
    BuildClassifier b = makeClassifier();

    Dataset data = TrainTestSplitTest.makeCircleDataset(new PortableRandom(42), 100, 200);
    TrainTestSplitTest.nukeData(data, 0.3);
    b.build(data);

    PredictClassifier p = b.getClassifier();

    //System.err.println(p.toString(new StringBuilder(), ""));

    SimpleEvaluation eval = new SimpleEvaluation();
    eval.evaluate(p, data);
    assertEquals(0.8433, eval.accuracy(), 0.01);
  }

  public void testSplitPointDoubleAveraging() {
    final Double prev = 0.38999999999999996;
    assertEquals(prev, RandomTreeBuilder.getSplitPoint(MlDataType.DOUBLE, prev, 0.39));
  }
}
