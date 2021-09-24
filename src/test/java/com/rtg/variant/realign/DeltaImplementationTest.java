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

import com.rtg.variant.util.arithmetic.LogPossibility;

/**
 */
public class DeltaImplementationTest extends AbstractRealignTest {

  @Override
  protected Delta getRealign(final Environment env, final RealignParams params) {
    final DeltaImplementation realign = new DeltaImplementation(LogPossibility.SINGLETON, params);
    realign.setEnv(env);
    realign.globalIntegrity();
    return realign;
  }

  /**
   * Test method for {@link com.rtg.util.integrity.IntegralAbstract#toString()}.
   */
  //  public void testToString() {
  //    TestUtils.containsAll(mRealign.toString(), "RealignImplementation(maxShift=2)", "ScoreMatrix");
  //  }

  public void testCombined() {
    final EnvironmentImplementation env = getEnv(TEMPLATE, READ);
    final Delta realign2 = getRealign(env);
    final DeltaImplementation realign = (DeltaImplementation) realign2;
    final String exp = ""
        + "Combined                   |G            0       |G            1       |A            2       |C            3       |G            4       |T            5       |A            6       |C            7       |G            8       |T            9       |T           10       |T           11       |N           12       |" + LS
        + "[  0]          0.001       |                     |         0.998  0.001|                     |                     |" + LS
        + "[  1]A                     |         0.001       |                     |         0.999       |                     |                     |" + LS
        + "[  2]C                     |                     |         0.001       |                0.002|         0.997       |                     |                     |" + LS
        + "[  3]G                     |                     |                     |         0.001  0.002|                0.005|         0.992       |                     |                     |" + LS
        + "[  4]C                     |                     |                     |                     |         0.003  0.002|         0.003  0.495|         0.498       |                     |                     |" + LS
        + "[  5]G                     |                     |                     |                     |                     |         0.004       |  0.004  0.498  0.495|         0.003       |                     |                     |" + LS
        + "[  6]A                     |                     |                     |                     |                     |                     |                     |         0.997  0.003|                     |                     |                     |" + LS
        + "[  7]C                     |                     |                     |                     |                     |                     |                     |                     |         1.000       |                     |                     |                     |" + LS
        + "[  8]G                     |                     |                     |                     |                     |                     |                     |                     |                     |         1.000       |                     |                     |                     |" + LS
        + "                           |G            0       |G            1       |A            2       |C            3       |G            4       |T            5       |A            6       |C            7       |G            8       |T            9       |T           10       |T           11       |N           12       |" + LS
        ;
    assertEquals(exp, realign.combinedToString());
  }

  public void testCombinedTerse() {
    final EnvironmentImplementation env = getEnv(TEMPLATE, READ);
    final Delta realign = getRealign(env);
    final String exp = ""
        + "         G   G   A   C   G   T   A   C   G   T   T   T   N   " + LS
        + "[  0]  . |   | 9.|   |   |" + LS
        + "[  1]A   | . |   | 9 |   |   |" + LS
        + "[  2]C   |   | . |  .| 9 |   |   |" + LS
        + "[  3]G   |   |   | ..|  .| 9 |   |   |" + LS
        + "[  4]C   |   |   |   | ..| .5| 5 |   |   |" + LS
        + "[  5]G   |   |   |   |   | . |.55| . |   |   |" + LS
        + "[  6]A   |   |   |   |   |   |   | 9.|   |   |   |" + LS
        + "[  7]C   |   |   |   |   |   |   |   | * |   |   |   |" + LS
        + "[  8]G   |   |   |   |   |   |   |   |   | * |   |   |   |" + LS
        + "         G   G   A   C   G   T   A   C   G   T   T   T   N   " + LS
        ;
    assertEquals(exp, realign.combinedTerse(env));
  }
}
