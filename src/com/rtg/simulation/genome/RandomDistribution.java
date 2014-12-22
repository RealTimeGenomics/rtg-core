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
package com.rtg.simulation.genome;

import com.rtg.util.PortableRandom;

/**
 * Produces random values with specified relative frequency values will be
 * between 0 (inclusive) and the size number of frequencies provided
 * (exclusive).
 *
 *
 */
public class RandomDistribution {

  private final int[] mDistribution;
  private final int mTotal;
  private final PortableRandom mGenerator;

  /**
   *
   * @param distribution relative distribution of the values (0, array length -
   *        1]
   * @param generator random number generator
   */
  public RandomDistribution(int[] distribution, PortableRandom generator) {
    mGenerator = generator;
    if (distribution != null) {
      mDistribution =  distribution.clone();
      mTotal = arraySum(distribution);
    } else {
      mDistribution = null;
      mTotal = 0;
    }
  }

  private int arraySum(int[] source) {
    int sum = 0;
    for (final int val : source) {
      sum += val;
    }
    return sum;
  }

  /**
   * @return a the next random value
   */
  public int nextValue() {
    if (mDistribution == null) {
      return 0;
    }
    int random = mGenerator.nextInt(mTotal);
    int i = 0;
    while (i < mDistribution.length && random >= 0) {
      random -= mDistribution[i];
      i++;
    }
    return i - 1;
  }

  /**
   * get the number of values in the distribution
   *
   * @return number of values in the distribution
   */
  public int valueCount() {
    return mDistribution.length;
  }

}
