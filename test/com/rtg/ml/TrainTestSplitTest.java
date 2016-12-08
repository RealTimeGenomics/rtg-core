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

import java.util.HashSet;

import com.rtg.util.PortableRandom;

import junit.framework.TestCase;

/**
 */
public class TrainTestSplitTest extends TestCase {


  static Dataset makeSimpleDataset(int posSize, int negSize) {
    final Dataset d = new Dataset(new Attribute("att1", MlDataType.DOUBLE));
    for (int i = 0; i < posSize; ++i) {
      d.addInstance(new Instance(new double[]{(double) i}, true));
    }
    for (int i = 0; i < negSize ; ++i) {
      d.addInstance(new Instance(new double[]{(double) (i + posSize)}, false));
    }
    return d;
  }

  // Positive examples are all inside a circle with radius 1
  static Dataset makeCircleDataset(PortableRandom random, int posSize, int negSize) {
    final Dataset d = new Dataset(new Attribute("X", MlDataType.DOUBLE), new Attribute("Y", MlDataType.DOUBLE));
    //double maxAngle = Math.PI / 2;
    final double maxAngle = Math.PI * 2;
    for (int i = 0; i < posSize; ++i) {
      final double angle = random.nextDouble() * maxAngle;
      final double radius = random.nextDouble();
      d.addInstance(new Instance(new double[]{radius * Math.cos(angle), radius * Math.sin(angle)}, true));
    }
    for (int i = 0; i < negSize ; ++i) {
      final double angle = random.nextDouble() * maxAngle;
      final double radius = 1.0 + random.nextDouble();
      d.addInstance(new Instance(new double[]{radius * Math.cos(angle), radius * Math.sin(angle)}, false));
    }
    return d;
  }

  public static void nukeData(Dataset data, double errorRate) {
    final PortableRandom rand = new PortableRandom(42);
    for (Instance inst : data.getInstances()) {
      for (int i = 0; i < inst.instance().length; ++i) {
        if (rand.nextDouble() < errorRate) {
          inst.instance()[i] = Double.NaN;
        }
      }
    }
  }

  public void runCombo(int numPos, int numNeg, int seed) {
    final Dataset d = makeSimpleDataset(numPos, numNeg);

    final int tsize = d.size() * 2 / 3;
    final TrainTestSplit t = TrainTestSplit.sampleWithReplacement(d, tsize, new PortableRandom(seed));
    assertEquals(tsize, t.mTrain.size());
    assertTrue(t.mTest.size() >= d.size() - tsize);

    final HashSet<Instance> train = new HashSet<>();
    train.addAll(t.mTrain.getInstances());

    assertTrue(train.size() > 0); // Chance of this failing at random is very low.

    for (Instance inst : t.mTest.getInstances()) {
      assertFalse(train.contains(inst));
    }
  }

  public void testNoddy() {
    runCombo(2, 1, 42);
  }

  public void testNoddy2() {
    try {
      runCombo(0, 0, 42);
      fail();
    } catch (IllegalArgumentException e) {
      // espected
    }
  }

  public void testMore() {
    runCombo(100, 100, 42);
  }
}
