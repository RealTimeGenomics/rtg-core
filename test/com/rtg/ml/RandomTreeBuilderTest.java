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

import java.util.Properties;

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
    final BuildClassifier b = makeClassifier();
    final Dataset data = TrainTestSplitTest.makeSimpleDataset(200, 100);

    final TrainTestSplit split = TrainTestSplit.sampleWithReplacement(data, 150, new PortableRandom(42));
    b.build(split.mTrain);

    final PredictClassifier p = b.getClassifier();

    final String thetree = p.toString(new StringBuilder(), "", data).toString();
    assertNotNull(thetree);
    //System.err.println(thetree);
  }

  public void testCircleTree() {
    final Dataset data = TrainTestSplitTest.makeCircleDataset(new PortableRandom(42), 300, 600);
    final Dataset testdata = TrainTestSplitTest.makeCircleDataset(new PortableRandom(92), 300, 600);
    TrainTestSplitTest.nukeData(testdata, 0.35, 0.15, Double.NaN);
    //TrainTestSplitTest.nukeData(testdata, 0.35, 0.15, 10);
    final Properties props = new Properties();

    final BuildClassifier b = makeClassifier();
    b.setProperties(props);
    buildAndEval(b, data, testdata, 0.89);
    //assertEquals(0.987, eval.accuracy(), 0.01);

    TrainTestSplitTest.nukeData(data, 0.35, 0.15, Double.NaN);
    config(b, "false", "false", "false");
    buildAndEval(b, data, testdata, 0.87);
    //assertEquals(0.8433, eval.accuracy(), 0.01);

    config(b, "true", "false", "false");
    buildAndEval(b, data, testdata, 0.84);
    //assertEquals(0.8433, eval.accuracy(), 0.01);

    config(b, "true", "false", "true");
    buildAndEval(b, data, testdata, 0.84);
    //assertEquals(0.8433, eval.accuracy(), 0.01);

    config(b, "true", "true", "false");
    buildAndEval(b, data, testdata, 0.79);
    //assertEquals(0.8433, eval.accuracy(), 0.01);

    config(b, "true", "true", "true");
    buildAndEval(b, data, testdata, 0.85);
    //assertEquals(0.8433, eval.accuracy(), 0.01);

  }

  private void config(BuildClassifier b, String ent, String prop, String split) {
    final Properties props = new Properties();
    props.setProperty(RandomTreeBuilder.PROP_ENTROPY_MISSING, ent);
    props.setProperty(RandomTreeBuilder.PROP_PROPAGATE_MISSING, prop);
    props.setProperty(RandomTreeBuilder.PROP_SPLIT_MISSING, split);
    b.setProperties(props);
  }

  private void buildAndEval(BuildClassifier builder, Dataset traindata, Dataset testdata, double expect) {
    builder.build(traindata);
    final PredictClassifier p = builder.getClassifier();
    final SimpleEvaluation eval = new SimpleEvaluation();
    eval.evaluate(p, testdata);
    //System.err.println("\n" + builder.toString());
    //System.err.print(p.toString(new StringBuilder(), "", data));
    //System.err.println("Nodes in tree: " + TestUtils.splitLines(p.toString(new StringBuilder(), "", traindata).toString()).length);
    //System.err.println("Accuracy: " + eval.accuracy());
    assertEquals(expect, eval.accuracy(), 0.01);
  }

  public void testSplitPointDoubleAveraging() {
    final Double prev = 0.38999999999999996;
    assertEquals(prev, RandomTreeBuilder.getSplitPoint(MlDataType.DOUBLE, prev, 0.39));
  }
}
