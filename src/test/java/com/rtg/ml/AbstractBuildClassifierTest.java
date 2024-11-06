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

import com.rtg.AbstractTest;
import com.rtg.util.PortableRandom;
import com.rtg.util.ThreadAware;

/**
 */
public abstract class AbstractBuildClassifierTest extends AbstractTest {

  abstract BuildClassifier makeClassifier();


  public void test() {
    final Attribute a = new Attribute("higgs", MlDataType.DOUBLE);
    final Dataset d = new Dataset(a);
    d.addInstance(new Instance(new double[] {Math.PI}, true));
    d.addInstance(new Instance(new double[] {0.0}, true));
    d.addInstance(new Instance(new double[] {4.2}, false));

    final BuildClassifier b = makeClassifier();
    b.build(d);

    final PredictClassifier p = b.getClassifier();

    for (Instance inst : d.getInstances()) {
      final double prob = p.predict(inst.instance());
      assertTrue(prob >= 0 && prob <= 1.0);
    }
  }

  public void testProperties() {
    final BuildClassifier b = makeClassifier();
    // Check that empty properties do not blow things up
    b.setProperties(new Properties());
    final Dataset data = TrainTestSplitTest.makeSimpleDataset(100, 100);
    b.build(data);
  }

  public void testMultithread() {
    final BuildClassifier b = makeClassifier();
    final Dataset data = TrainTestSplitTest.makeSimpleDataset(100, 100);

    final TrainTestSplit split = TrainTestSplit.sampleWithReplacement(data, 150, new PortableRandom(42));
    b.build(split.mTrain);

    final PredictClassifier p = b.getClassifier();

    //System.err.println(p.toString(new StringBuilder(), ""));

    assertTrue(p == b.getClassifier()); // You can ask for the same classifier multiple times

    if (b instanceof ThreadAware) {
      for (int threads = 1; threads < 10; threads += 3) {
        ((ThreadAware) b).setNumberOfThreads(threads);
        b.build(split.mTrain);
        final PredictClassifier p2 = b.getClassifier();

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
