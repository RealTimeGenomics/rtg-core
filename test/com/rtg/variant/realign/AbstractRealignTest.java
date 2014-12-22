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

package com.rtg.variant.realign;

import static com.rtg.util.StringUtils.LS;

import com.rtg.mode.DnaUtils;
import com.rtg.util.StringUtils;
import com.rtg.variant.realign.ScoreMatrixTest.MockRealignParams;

import junit.framework.TestCase;


/**
 */
public abstract class AbstractRealignTest extends TestCase {

  static final byte[] READ = DnaUtils.encodeString(ScoreMatrixTest.READ);
  static final byte[] TEMPLATE = DnaUtils.encodeString(ScoreMatrixTest.TEMPLATE);

  protected abstract Delta getRealign(final Environment env, final RealignParams params);

  protected Delta getRealign(final EnvironmentImplementation env) {
    final RealignParams params = new MockRealignParams();
    return getRealign(env, params);
  }

  protected EnvironmentImplementation getEnv(final byte[] template, final byte[] read) {
    final double[] qualities = new double[read.length];
    qualities[0] = 0.1;
    qualities[1] = Math.pow(10.0, -15 / 10.0);
    for (int i = 2; i < qualities.length; i++) {
      qualities[i] = 0.01;
    }
    final EnvironmentImplementation env = new EnvironmentImplementation(
        2, //maxShift
        template,
        2, //start
        read,
        qualities
        );
    env.integrity();
    return env;
  }

  public void testRealign() {
    final EnvironmentImplementation env = getEnv(TEMPLATE, READ);
    final Delta realign = getRealign(env);
    assertNotNull(realign);
    realign.globalIntegrity();
  }

  /**
   * Test method for {@link com.rtg.variant.realign.DeltaImplementation#totalScoreLn()}.
   */
  public void testLnScore() {
    final EnvironmentImplementation env = getEnv(TEMPLATE, READ);
    final Delta realign = getRealign(env);
    assertEquals(-13.2824, realign.totalScoreLn(), 0.0002);
  }

  /**
   * The expected probabilities are calculated in scorematrix.xls.
   */
  public void testProbabilitiesLn() {
    final EnvironmentImplementation env = getEnv(TEMPLATE, READ);
    final Delta realign = getRealign(env);
    double[] probs;
    // the expected start of the read, but this is before the best
    // alignment, so no strong evidence.
    probs = realign.probabilitiesLn(2);
    assertNotNull(probs);
    assertEquals(4, probs.length);
    assertEquals(-13.2824, probs[0], 0.0001);
    assertEquals(-16.3197, probs[1], 0.0001);
    assertEquals(-15.1021, probs[2], 0.0001);
    assertEquals(-16.4591, probs[3], 0.0001);

    // an ambiguous position, where C and G are equally likely
    probs = realign.probabilitiesLn(5);
    assertNotNull(probs);
    assertEquals(-13.2761, probs[0], 0.0001);
    assertEquals(-8.9356, probs[1], 0.0001);
    assertEquals(-8.9353, probs[2], 0.0001);
    assertEquals(-13.2824, probs[3], 0.0001);

    // just past the end of the best alignment
    probs = realign.probabilitiesLn(9);
    assertNotNull(probs);
    assertEquals(-13.2801, probs[0], 0.0001);
    assertEquals(-13.2824, probs[1], 0.0001);
    assertEquals(-13.2809, probs[2], 0.0001);
    assertEquals(-13.2825, probs[3], 0.0001);
  }

  public void testWithinRead() {
    final EnvironmentImplementation env = getEnv(TEMPLATE, READ);
    final Delta realign = getRealign(env);
    // expected start is 2.
    assertEquals(Double.NEGATIVE_INFINITY, realign.withinReadLn(-1));
    assertEquals(-7.17996, realign.withinReadLn(0), 0.0001);
    assertEquals(-7.03252, realign.withinReadLn(1), 0.0001);
    assertEquals(-0.00053, realign.withinReadLn(2), 0.0001);
    assertEquals(-0.00016, realign.withinReadLn(3), 0.0001);
    assertEquals(0.0, realign.withinReadLn(4));
    assertEquals(0.0, realign.withinReadLn(5));
    assertEquals(0.0, realign.withinReadLn(6));
    assertEquals(0.0, realign.withinReadLn(7));
    assertEquals(-0.0002, realign.withinReadLn(8), 0.0001);
    // expected end is 2 + (8 - 1) = 9, but actual end is 8.
    assertEquals(-10.5400, realign.withinReadLn(9), 0.0001);
    assertEquals(-11.2175, realign.withinReadLn(10), 0.0001);
    assertEquals(-11.3468, realign.withinReadLn(11), 0.0001);
    assertEquals(Double.NEGATIVE_INFINITY, realign.withinReadLn(12));
  }

  void checkPredictions(final String template, final String read, final String expected) {
    final StringBuilder sb = new StringBuilder();
    final EnvironmentImplementation env = getEnv(DnaUtils.encodeString(template), DnaUtils.encodeString(read));
    final Delta realign = getRealign(env);
    for (int i = 0; i < template.length(); i++) {
      final ProbabilityArray pr = new ProbabilityArray(realign.probabilitiesLn(i));
      sb.append("[").append(i).append("]");
      sb.append(pr);
      sb.append(StringUtils.LS);

    }
    assertEquals(expected, sb.toString());
  }

  public void testSnps0() {
    //the case from the spreadsheet - numbers checked against the spreedsheet
    final String exp = ""
        + "[0][0.253 0.249 0.249 0.249]" + LS
        + "[1][0.247 0.260 0.247 0.247]" + LS
        + "[2][0.799 0.038 0.129 0.033]" + LS
        + "[3][0.014 0.945 0.027 0.014]" + LS
        + "[4][0.009 0.337 0.648 0.006]" + LS
        + "[5][0.006 0.494 0.494 0.006]" + LS
        + "[6][0.655 0.008 0.331 0.006]" + LS
        + "[7][0.010 0.976 0.007 0.006]" + LS
        + "[8][0.007 0.007 0.979 0.007]" + LS
        + "[9][0.250 0.250 0.250 0.250]" + LS
        + "[10][0.250 0.250 0.250 0.250]" + LS
        + "[11][0.250 0.250 0.250 0.250]" + LS
        ;
    checkPredictions(ScoreMatrixTest.TEMPLATE, ScoreMatrixTest.READ, exp);
    //checkPredictions("AAACGCGACGTT", "ACGCGACG", "");
  }

}
