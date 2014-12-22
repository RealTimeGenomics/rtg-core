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
    for (int i = 0; i < subsetSize; i++) {
      final int index = seed.nextInt(size);
      train.addInstance(input.getInstances().get(index));
      used[index] = true;
    }

    // Put any instances not used for training into split.mTest
    final Dataset test = new Dataset(input.getAttributes());
    for (int i = 0; i < size; i++) {
      if (!used[i]) {
        test.addInstance(input.getInstances().get(i));
      }
    }

    return new TrainTestSplit(train, test);
  }
}
