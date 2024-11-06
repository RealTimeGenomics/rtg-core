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

  public void runCombo(int numPos, int numNeg, int seed) {
    final Dataset d = makeSimpleDataset(numPos, numNeg);

    final int tsize = d.size() * 2 / 3;
    final TrainTestSplit t = TrainTestSplit.sampleWithReplacement(d, tsize, new PortableRandom(seed));
    assertEquals(tsize, t.mTrain.size());
    assertTrue(t.mTest.size() >= d.size() - tsize);

    final HashSet<Instance> train = new HashSet<>(t.mTrain.getInstances());

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
