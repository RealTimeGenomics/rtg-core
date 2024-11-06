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
