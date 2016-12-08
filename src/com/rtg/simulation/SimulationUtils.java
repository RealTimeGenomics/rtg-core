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
package com.rtg.simulation;

import java.util.Arrays;

/**
 * Utility methods for simulation
 */
public final class SimulationUtils {

  private SimulationUtils() { }

  /**
   * Generates an accumulated distribution from an input distribution. The input distribution need not sum to 1
   * @param dist the non-cumulative distribution
   * @return the cumulative distribution
   */
  public static double[] cumulativeDistribution(final double... dist) {
    double sum = 0;
    for (final double r : dist) {
      sum += r;
    }

    final double[] thres = new double[dist.length];
    double currentThres = 0;
    for (int i = 0; i < dist.length; ++i) {
      currentThres += dist[i];
      thres[i] = currentThres / sum;
    }
    return thres;
  }

  /**
   * Find entry position in a cumulative distribution.
   * @param dist cumulative distribution
   * @param rand a double chosen between 0.0 and 1.0
   * @return a chosen length
   */
  public static int chooseLength(final double[] dist, final double rand) {
    assert rand <= 1.0 && rand >= 0.0;
    int len = Arrays.binarySearch(dist, rand);
    if (len < 0) {
      len = -len - 1;
    }
    while (len < dist.length && dist[len] == 0.0) {
      ++len;
    }
    return len;
  }

}
