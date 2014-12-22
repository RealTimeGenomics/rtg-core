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
import com.rtg.util.ThreadAware;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractBuildClassifierTest extends TestCase {

  abstract BuildClassifier makeClassifier();


  public void test() {
    final Attribute a = new Attribute("higgs", MlDataType.DOUBLE);
    final Dataset d = new Dataset(a);
    d.addInstance(new Instance(new double[] {Math.PI}, true));
    d.addInstance(new Instance(new double[] {0.0}, true));
    d.addInstance(new Instance(new double[] {4.2}, false));

    BuildClassifier b = makeClassifier();
    b.build(d);

    PredictClassifier p = b.getClassifier();

    for (Instance inst : d.getInstances()) {
      double prob = p.predict(inst.instance());
      assertTrue(prob >= 0 && prob <= 1.0);
    }
  }

  public void testProperties() {
    BuildClassifier b = makeClassifier();
    // Check that empty properties do not blow things up
    b.setProperties(new Properties());
    Dataset data = TrainTestSplitTest.makeSimpleDataset(100, 100);
    b.build(data);
  }

  public void testMultithread() {
    BuildClassifier b = makeClassifier();
    Dataset data = TrainTestSplitTest.makeSimpleDataset(100, 100);

    TrainTestSplit split = TrainTestSplit.sampleWithReplacement(data, 150, new PortableRandom(42));
    b.build(split.mTrain);

    PredictClassifier p = b.getClassifier();

    //System.err.println(p.toString(new StringBuilder(), ""));

    assertTrue(p == b.getClassifier()); // You can ask for the same classifier multiple times

    if (b instanceof ThreadAware) {
      for (int threads = 1; threads < 10; threads += 3) {
        ((ThreadAware) b).setNumberOfThreads(threads);
        b.build(split.mTrain);
        PredictClassifier p2 = b.getClassifier();

        assertFalse(p == p2);             // But every build must create a new classifier

        // We should get the same predictions as the single-thread classifier
        for (Instance inst : split.mTest.getInstances()) {
          final double predict = p.predict(inst.instance());
          final double predict2 = p2.predict(inst.instance());
          assertEquals("Threads=" + threads + " " + predict + " != " + predict2, predict, predict2);
        }
      }
    }
  }

}
