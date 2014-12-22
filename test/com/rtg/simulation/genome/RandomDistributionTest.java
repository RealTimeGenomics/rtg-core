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

import junit.framework.TestCase;

/**
 *
 */
public class RandomDistributionTest extends TestCase {

  /**
   */
  public RandomDistributionTest(final String name) {
    super(name);
  }

  public void testFrequencies() {
    final PortableRandom rand = new PortableRandom(2);
    final int[] distribution = {1, 2, 6, 4, 0, 9};
    final RandomDistribution dist = new RandomDistribution(distribution, rand);
    final int[] results = new int[distribution.length];
    for (int i = 0; i < results.length; i++) {
      results[i] = 0;
    }
    for (int i = 0; i < 1000000; i++) {
      results[dist.nextValue()]++;
    }
    assertEquals(2, ((double) results[1]) / results[0], 0.1);
    assertEquals(6, ((double) results[2]) / results[0], 0.1);
    assertEquals(4, ((double) results[3]) / results[0], 0.1);
    assertEquals(0, results[4]);
    assertEquals(9, ((double) results[5]) / results[0], 0.1);

  }

  public void testZeroDistribution() {
    final PortableRandom rand = new PortableRandom(2);
    final int[] distribution = {0, 0, 0, 1};
    final RandomDistribution dist = new RandomDistribution(distribution, rand);
    final int[] results = new int[distribution.length];
    for (int i = 0; i < results.length; i++) {
      results[i] = 0;
    }
    for (int i = 0; i < 1000000; i++) {
      results[dist.nextValue()]++;
    }
    assertEquals(0, results[0]);
    assertEquals(0, results[1]);
    assertEquals(0, results[2]);
    assertEquals(1000000, results[3]);

  }

  public void testValueCount() {
    final PortableRandom rand = new PortableRandom(1);
    final int[] distribution = {1, 2, 6, 4, 9};
    final RandomDistribution dist = new RandomDistribution(distribution, rand);
    assertEquals(5, dist.valueCount());
  }

}
