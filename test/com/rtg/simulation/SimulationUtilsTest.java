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

import junit.framework.TestCase;

/**
 */
public class SimulationUtilsTest extends TestCase {



  public void testDistUtils() {

    final double[] cumDist = SimulationUtils.cumulativeDistribution(0.1, 0.1, 0.1, 0.4, 0.1);

    final double[] expected = {0.125, 0.25, 0.375, 0.875, 1.0};
    assertEquals(expected.length, cumDist.length);
    for (int i = 0; i < expected.length; ++i) {
      assertEquals(expected[i], cumDist[i], 0.000001);
    }

    assertEquals(0, SimulationUtils.chooseLength(cumDist, 0.1));
    assertEquals(1, SimulationUtils.chooseLength(cumDist, 0.25));
    assertEquals(2, SimulationUtils.chooseLength(cumDist, 0.3));
    assertEquals(3, SimulationUtils.chooseLength(cumDist, 0.8));
    assertEquals(4, SimulationUtils.chooseLength(cumDist, 1.0));
  }

}
