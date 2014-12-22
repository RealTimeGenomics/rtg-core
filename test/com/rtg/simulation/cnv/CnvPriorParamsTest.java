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
package com.rtg.simulation.cnv;

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.util.Arrays;

import com.rtg.simulation.SimulationUtils;
import com.rtg.simulation.cnv.CnvPriorParams.CnvPriorParamsBuilder;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class CnvPriorParamsTest extends TestCase {

  /**
   */
  public CnvPriorParamsTest(String name) {
    super(name);
  }

  private final class MyCnvPriorParams extends CnvPriorParams {

    /**
     * Create params
     * @param builder params builder
     * @throws InvalidParamsException
     */
    public MyCnvPriorParams(CnvPriorParamsBuilder builder) {
      super(builder);
      try {
        checkEndsList(new int[] {5, 1, 10});
        fail();
      } catch (IllegalArgumentException e) {
        TestUtils.containsAll(e.getMessage(), "ends list must contain values > 1 with from left to right increasing or equal values");
      }
      try {
        checkEndsList(new int[] {0, 1, 10});
        fail();
      } catch (IllegalArgumentException e) {
        TestUtils.containsAll(e.getMessage(), "ends list must contain values > 1 with from left to right increasing or equal values");
      }
     try {
        checkDistribution(new double[] {0.5, 0.5, 0.5, 0.6});
        fail();
      } catch (IllegalArgumentException e) {
        TestUtils.containsAll(e.getMessage(), "distribution must sum to 1.0, not 2");
      }
      try {
        checkProbability(1.1);
        fail();
      } catch (IllegalArgumentException e) {
        TestUtils.containsAll(e.getMessage(), "rate must be 0.0 .. 1.0, not 1.1");
      }
      double[] powerDist = powerLengthDistribution();
      double[] powerThres = SimulationUtils.cumulativeDistribution(powerDist);
      assertTrue(Arrays.equals(powerThres, powerLengthThresholds()));

      double[][] cnDist = copyNumberDistribution();
      double[][] cnThres = copyNumberThresholds();
      double[] thres = SimulationUtils.cumulativeDistribution(cnDist[0]);
      assertTrue(Arrays.equals(thres, cnThres[0]));
      thres = SimulationUtils.cumulativeDistribution(cnDist[1]);
      assertTrue(Arrays.equals(thres, cnThres[1]));
      thres = SimulationUtils.cumulativeDistribution(cnDist[2]);
      assertTrue(Arrays.equals(thres, cnThres[2]));
      thres = SimulationUtils.cumulativeDistribution(cnDist[3]);
      assertTrue(Arrays.equals(thres, cnThres[3]));
      thres = SimulationUtils.cumulativeDistribution(cnDist[4]);
      assertTrue(Arrays.equals(thres, cnThres[4]));
      thres = SimulationUtils.cumulativeDistribution(cnDist[5]);
      assertTrue(Arrays.equals(thres, cnThres[5]));
      thres = SimulationUtils.cumulativeDistribution(cnDist[6]);
      assertTrue(Arrays.equals(thres, cnThres[6]));
      thres = SimulationUtils.cumulativeDistribution(cnDist[7]);
      assertTrue(Arrays.equals(thres, cnThres[7]));
    }
  }

  private CnvPriorParams mPriors;

  @Override
  public void setUp() throws Exception {
    mPriors = CnvPriorParams.builder().cnvpriors("testcnv-default").create();
  }

  @Override
  public void tearDown() {
    mPriors = null;
  }

  public void testMyCnvPriorParams() {
    new MyCnvPriorParams(CnvPriorParams.builder());
  }

  public void test() throws InvalidParamsException, IOException {
    mPriors = CnvPriorParams.builder()
    .cnvpriors("testcnv-default")
    .create();
    assertEquals(0.225, mPriors.probDeletedOnOneStrand());
    assertEquals(0.075, mPriors.probDeletedOnBothStrands());
  }

  /** Check that the property values are read correctly. */
  public void testDefaults() {
    final double[] powerDist = mPriors.powerLengthDistribution();
    assertEquals(8, powerDist.length);
    assertEquals(0.0, powerDist[0]);
    assertEquals(0.0, powerDist[1]);
    assertEquals(0.5, powerDist[2]);
    assertEquals(0.2, powerDist[3]);
    assertEquals(0.1, powerDist[4]);
    assertEquals(0.1, powerDist[5]);
    assertEquals(0.05, powerDist[6]);
    assertEquals(0.05, powerDist[7]);

    final int[] endsList = mPriors.copyRangeEnds();
    assertTrue(Arrays.equals(endsList, new int[]{1, 2, 5, 10, 20, 100}));


    final double[][] cnDist = mPriors.copyNumberDistribution();
    assertTrue(Arrays.equals(cnDist[0], new double[] {0.4, 0.3, 0.1, 0.2, 0.0, 0.0}));
    assertTrue(Arrays.equals(cnDist[1], new double[] {0.4, 0.3, 0.1, 0.2, 0.0, 0.0}));
    assertTrue(Arrays.equals(cnDist[2], new double[] {0.4, 0.3, 0.1, 0.2, 0.0, 0.0}));
    assertTrue(Arrays.equals(cnDist[3], new double[] {0.4, 0.3, 0.1, 0.2, 0.0, 0.0}));
    assertTrue(Arrays.equals(cnDist[4], new double[] {0.6, 0.3, 0.1, 0.0, 0.0, 0.0}));
    assertTrue(Arrays.equals(cnDist[5], new double[] {0.6, 0.3, 0.1, 0.0, 0.0, 0.0}));
    assertTrue(Arrays.equals(cnDist[6], new double[] {0.6, 0.3, 0.1, 0.0, 0.0, 0.0}));
    assertTrue(Arrays.equals(cnDist[7], new double[] {0.6, 0.3, 0.1, 0.0, 0.0, 0.0}));
  }

  public void testBuilder() {
    final CnvPriorParamsBuilder builder = CnvPriorParams.builder();
    assertNotNull(builder);
    assertTrue(builder == builder.copyRangeEnds(new int[] {1, 3}));
  }

  public void testToString() {
    assertEquals("    probability one delete = 0.225" + LS
        + "probability both delete = 0.075" + LS,
        mPriors.toString());
  }
}
