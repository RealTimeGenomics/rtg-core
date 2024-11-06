/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    testdata.injectErrors(0.35, 0.15, Double.NaN);
    //TrainTestSplitTest.nukeData(testdata, 0.35, 0.15, 10);
    final Properties props = new Properties();

    final BuildClassifier b = makeClassifier();
    b.setProperties(props);
    buildAndEval(b, data, testdata, 0.89);

    data.injectErrors(0.35, 0.15, Double.NaN);
    config(b, "false", "false", "false");
    buildAndEval(b, data, testdata, 0.87);

    config(b, "true", "false", "false");
    buildAndEval(b, data, testdata, 0.84);

    config(b, "true", "false", "true");
    buildAndEval(b, data, testdata, 0.84);

    config(b, "true", "both", "false");
    buildAndEval(b, data, testdata, 0.79);

    config(b, "true", "both", "true");
    buildAndEval(b, data, testdata, 0.85);

    config(b, "true", "random", "false");
    buildAndEval(b, data, testdata, 0.81);

    config(b, "true", "random", "true");
    buildAndEval(b, data, testdata, 0.85);
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
