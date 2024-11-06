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

package com.rtg.variant.bayes;

import static com.rtg.util.StringUtils.LS;

import com.rtg.reference.Ploidy;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.HypothesesMock;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class HypothesesTest extends TestCase {

  private void checkName(final Hypotheses<?> hy, final String exp, final int n) {
    assertEquals(exp, hy.name(n));
  }

  public void testDiploid() {
    final DescriptionCommon description = DescriptionSnp.SINGLETON;
    final PossibilityArithmetic arithmetic = SimplePossibility.SINGLETON;
    final Hypotheses<?> hy = new MockHypotheses<>(description, arithmetic, false, null, 0);
    assertTrue(description == hy.description());
    assertTrue(arithmetic == hy.arithmetic());
    assertEquals(10, hy.code().size());
    assertEquals(10, hy.size());
    assertFalse(hy.haploid());
    assertEquals(Ploidy.DIPLOID, hy.ploidy());

    assertTrue(hy.valid(0));
    assertTrue(hy.valid(9));
    assertFalse(hy.valid(-1));
    assertFalse(hy.valid(10));

    checkName(hy, "A:A", 0);
    checkName(hy, "T:T", 3);
    checkName(hy, "A:C", 4);
    checkName(hy, "C:G", 5);
    checkName(hy, "A:T", 9);

    assertTrue(hy.homozygous(0));
    assertTrue(hy.homozygous(3));
    assertFalse(hy.homozygous(4));
    assertFalse(hy.homozygous(9));
  }

  public void testHaploid() {
    final DescriptionCommon description = DescriptionSnp.SINGLETON;
    final PossibilityArithmetic arithmetic = SimplePossibility.SINGLETON;
    final Hypotheses<?> hy = new MockHypotheses<>(description, arithmetic, true, null, 0);
    assertTrue(description == hy.description());
    assertTrue(arithmetic == hy.arithmetic());
    assertEquals(4, hy.code().size());
    assertEquals(4, hy.size());
    assertTrue(hy.haploid());
    assertEquals(Ploidy.HAPLOID, hy.ploidy());

    assertTrue(hy.valid(0));
    assertTrue(hy.valid(3));
    assertFalse(hy.valid(-1));
    assertFalse(hy.valid(4));

    checkName(hy, "A", 0);
    checkName(hy, "T", 3);
  }

  public void testDifferentLengths() {
    final DescriptionCommon description = new DescriptionCommon("X", "YY");
    final PossibilityArithmetic arithmetic = SimplePossibility.SINGLETON;
    final Hypotheses<?> hy = new MockHypotheses<>(description, arithmetic, true, null, 0);
    assertTrue(description == hy.description());
    assertTrue(arithmetic == hy.arithmetic());
    assertEquals(2, hy.code().size());
    assertEquals(2, hy.size());
    assertTrue(hy.haploid());
  }

  public void testToStringDiploid() {
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final HypothesesPrior<Description> hy = new HypothesesMock<Description>(new DescriptionCommon("X", "YZ"), arith, false, 1) {
      @Override
      protected void initPriors(double[] priors) {
        for (int i = 0; i < priors.length; ++i) {
          priors[i] = arith.prob2Poss(1.0 / (i + 2));
        }
      }
    };
    assertEquals(5, hy.maxNameLength());
    final String exp = ""
        + "Hypotheses" + LS
        + "   X:X  -0.693" + LS
        + " YZ:YZ  -1.099" + LS
        + "  X:YZ  -1.386" + LS
        ;
    assertEquals(exp, hy.toString());
    assertEquals(0.5, hy.p(0));
    assertEquals(1.0 / 3.0, hy.p(1));
    assertEquals(0.25, hy.p(2));

  }

  public void testToStringHaploid() {
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Hypotheses<Description> hy = new HypothesesMock<Description>(new DescriptionCommon("X", "YZ"), arith, true, 1) {
      @Override
      protected void initPriors(double[] priors) {
        for (int i = 0; i < priors.length; ++i) {
          priors[i] = arith.prob2Poss(1.0 / (i + 2));
        }
      }
    };
    assertEquals(2, hy.maxNameLength());
    final String exp = ""
        + "Hypotheses" + LS
        + "  X  -0.693" + LS
        + " YZ  -1.099" + LS
        ;
    assertEquals(exp, hy.toString());

  }
}
