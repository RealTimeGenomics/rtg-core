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
import com.rtg.variant.realign.AbstractScoreMatrixCGTest.MockRealignParamsCG;
import com.rtg.variant.util.arithmetic.LogPossibility;

import junit.framework.TestCase;

/**
 */
public class CgRealignTest extends TestCase {

  /** This is the read and template from the <code>scorematrixCG.xls</code> spreadsheet */
  private static final byte[] READ = DnaUtils.encodeString(AbstractScoreMatrixCGTest.READ);
  private static final byte[] TEMPLATE = DnaUtils.encodeString(AbstractScoreMatrixCGTest.TEMPLATE);

  protected Delta getRealign(final Environment env, final RealignParams params) {
    final DeltaSlowly realign = new DeltaSlowly(
        LogPossibility.SINGLETON,
        new ScoreMatrixCG(LogPossibility.SINGLETON, params),
        new ScoreMatrixCG(LogPossibility.SINGLETON, params), new ScoreMatrixCGReverse(LogPossibility.SINGLETON, params)
        );
    realign.setEnv(env);
    //realign.globalIntegrity(); //TODO check this - see comment in globalIntegrity that it wont work for CG
    return realign;
  }

  protected Delta getRealign(final byte[] template, final byte[] read) {
    final double[] qualities = new double[read.length];
    qualities[0] = 0.1;
    qualities[1] = Math.pow(10.0, -15 / 10.0);
    for (int i = 2; i < qualities.length; ++i) {
      qualities[i] = 0.01;
    }
    final EnvironmentImplementation env = new EnvironmentImplementation(
        7, //maxShift
        template,
        7, //start
        read,
        qualities
        );
    env.integrity();
    final RealignParams params = new MockRealignParamsCG();
    return getRealign(env, params);
  }

  public void testRealign() {
    final Delta realign = getRealign(TEMPLATE, READ);
    assertNotNull(realign);
    // realign.globalIntegrity();  // TODO: make this work for CG
  }

  /**
   * Test method for {@link com.rtg.variant.realign.DeltaImplementation#totalScoreLn()}.
   */
  public void testLnScore() {
    final Delta realign = getRealign(TEMPLATE, READ);
    assertEquals(-10.01247, realign.totalScoreLn(), 0.0002);
  }

  public void testWithinRead() {
    final Delta realign = getRealign(TEMPLATE, READ);
    // expected start is 7.
    assertEquals(Double.NEGATIVE_INFINITY, realign.withinReadLn(-1));

    assertEquals(-23.9712, realign.withinReadLn(0), 0.0001);

    assertEquals(-14.0367, realign.withinReadLn(6), 0.0001);
    assertEquals(-4.5470, realign.withinReadLn(7), 0.0001);
    assertEquals(-4.5435, realign.withinReadLn(8), 0.0001);
    assertEquals(-2.3226, realign.withinReadLn(9), 0.0001);
    assertEquals(-0.0140, realign.withinReadLn(10), 0.0001);

    assertEquals(0.0, realign.withinReadLn(15));
    assertEquals(0.0, realign.withinReadLn(48), 0.0001);

    assertEquals(-0.0002, realign.withinReadLn(49), 0.0001);
    assertEquals(-7.6540, realign.withinReadLn(50), 0.0001);
    assertEquals(-14.1332, realign.withinReadLn(51), 0.0001);
    assertEquals(-15.9250, realign.withinReadLn(52), 0.0001);
    assertEquals(Double.NEGATIVE_INFINITY, realign.withinReadLn(53));
  }

  void checkPredictions(final String template, final String read, final String expected) {
    final StringBuilder sb = new StringBuilder();
    final Delta realign = getRealign(DnaUtils.encodeString(template), DnaUtils.encodeString(read));
    for (int i = 0; i < template.length(); ++i) {
      final ProbabilityArray pr = new ProbabilityArray(realign.probabilitiesLn(i));
      sb.append("[").append(i).append("]");
      sb.append(pr);
      sb.append(LS);

    }
    assertEquals(expected, sb.toString());
  }

  public void testSnps0() {
    // the case from the spreadsheet - these numbers look great,
    // but the spreadsheet doesn't calculate them yet.
    final String exp = ""
        + "[0][0.250 0.250 0.250 0.250]" + LS
        + "[1][0.250 0.250 0.250 0.250]" + LS
        + "[2][0.250 0.250 0.250 0.250]" + LS
        + "[3][0.250 0.250 0.250 0.250]" + LS
        + "[4][0.250 0.250 0.250 0.250]" + LS
        + "[5][0.250 0.250 0.250 0.250]" + LS
        + "[6][0.250 0.250 0.250 0.250]" + LS
        + "[7][0.252 0.249 0.249 0.249]" + LS
        + "[8][0.250 0.249 0.249 0.252]" + LS
        + "[9][0.269 0.244 0.244 0.244]" + LS
        + "[10][0.133 0.019 0.019 0.830]" + LS
        + "[11][0.015 0.014 0.014 0.958]" + LS
        + "[12][0.892 0.014 0.015 0.080]" + LS
        + "[13][0.989 0.004 0.004 0.004]" + LS
        + "[14][0.991 0.003 0.004 0.003]" + LS
        + "[15][0.007 0.006 0.980 0.006]" + LS
        + "[16][0.006 0.006 0.981 0.006]" + LS
        + "[17][0.006 0.981 0.006 0.006]" + LS
        + "[18][0.006 0.006 0.981 0.006]" + LS
        + "[19][0.981 0.006 0.006 0.006]" + LS
        + "[20][0.006 0.981 0.006 0.006]" + LS
        + "[21][0.980 0.006 0.006 0.007]" + LS
        + "[22][0.007 0.007 0.007 0.979]" + LS
        + "[23][0.007 0.007 0.980 0.007]" + LS
        + "[24][0.006 0.980 0.007 0.006]" + LS
        + "[25][0.006 0.981 0.006 0.006]" + LS
        + "[26][0.981 0.006 0.006 0.006]" + LS
        + "[27][0.981 0.006 0.006 0.006]" + LS
        + "[28][0.006 0.006 0.006 0.981]" + LS
        + "[29][0.006 0.006 0.981 0.006]" + LS
        + "[30][0.006 0.006 0.006 0.981]" + LS
        + "[31][0.006 0.006 0.981 0.006]" + LS
        + "[32][0.006 0.006 0.006 0.981]" + LS
        + "[33][0.250 0.250 0.250 0.250]" + LS
        + "[34][0.250 0.250 0.250 0.251]" + LS
        + "[35][0.250 0.250 0.250 0.250]" + LS
        + "[36][0.250 0.250 0.250 0.250]" + LS
        + "[37][0.250 0.250 0.250 0.250]" + LS
        + "[38][0.250 0.250 0.250 0.251]" + LS
        + "[39][0.250 0.250 0.250 0.251]" + LS
        + "[40][0.013 0.015 0.013 0.960]" + LS
        + "[41][0.010 0.011 0.010 0.969]" + LS
        + "[42][0.006 0.981 0.006 0.006]" + LS
        + "[43][0.981 0.006 0.006 0.006]" + LS
        + "[44][0.981 0.006 0.006 0.006]" + LS
        + "[45][0.006 0.981 0.006 0.006]" + LS
        + "[46][0.006 0.006 0.006 0.981]" + LS
        + "[47][0.006 0.006 0.006 0.981]" + LS
        + "[48][0.006 0.007 0.006 0.980]" + LS
        + "[49][0.007 0.978 0.007 0.008]" + LS
        + "[50][0.250 0.250 0.250 0.250]" + LS
        + "[51][0.250 0.250 0.250 0.250]" + LS
        + "[52][0.250 0.250 0.250 0.250]" + LS
        + "[53][0.250 0.250 0.250 0.250]" + LS
        + "[54][0.250 0.250 0.250 0.250]" + LS
        + "[55][0.250 0.250 0.250 0.250]" + LS
        + "[56][0.250 0.250 0.250 0.250]" + LS
        ;
    checkPredictions(AbstractScoreMatrixCGTest.TEMPLATE, AbstractScoreMatrixCGTest.READ, exp);
  }

}
