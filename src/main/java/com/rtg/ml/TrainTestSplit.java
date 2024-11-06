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

import com.rtg.util.PortableRandom;

/**
 * Encapsulates a training test set split.
 *
 */
public class TrainTestSplit {

  final Dataset mTrain;
  final Dataset mTest;

  TrainTestSplit(Dataset train, Dataset test) {
    mTrain = train;
    mTest = test;
  }

  /**
   * Create a training test split using replacement. The test set consists of all instances that are not used for training.
   * @param input the input dataset
   * @param subsetSize the target number of training instances
   * @param seed random seed
   * @return the TrainTestSplit
   */
  static TrainTestSplit sampleWithReplacement(Dataset input, int subsetSize, PortableRandom seed) {
    if (input.size() == 0) {
      throw new IllegalArgumentException("Dataset size is 0");
    }
    if (subsetSize <= 0 || subsetSize > input.size()) {
      throw new IllegalArgumentException("Subset size must be between 0 and " + input.size());
    }

    final int size = input.size();
    final boolean[] used = new boolean[size];

    // Choose subsetSize instances from d, poke into split.mTrain
    final Dataset train = new Dataset(input.getAttributes());
    for (int i = 0; i < subsetSize; ++i) {
      final int index = seed.nextInt(size);
      train.addInstance(input.getInstances().get(index));
      used[index] = true;
    }

    // Put any instances not used for training into split.mTest
    final Dataset test = new Dataset(input.getAttributes());
    for (int i = 0; i < size; ++i) {
      if (!used[i]) {
        test.addInstance(input.getInstances().get(i));
      }
    }

    return new TrainTestSplit(train, test);
  }
}
