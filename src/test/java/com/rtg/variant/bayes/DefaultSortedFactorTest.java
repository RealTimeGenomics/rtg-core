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

import com.rtg.variant.bayes.multisample.forwardbackward.MutableFactor;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class DefaultSortedFactorTest extends TestCase {

  public void test() {
    final Description desc = new DescriptionCommon("A", "B");
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final Hypotheses<Description> hyp = new MockHypotheses<>(desc, arith, false, new double[3], 0);
    final Factor<Description> m = new MutableFactor<>(hyp, arith, new double[] {0.4, 0.6, 1.0});
    assertEquals(3, m.size());
    assertEquals(0.4, m.p(0), 1e-10);
    assertEquals(0.6, m.p(1), 1e-10);
    assertEquals(1.0, m.p(2), 1e-10);
    assertFalse(m.isNormalized());
    final SortedFactor<Description> s = new DefaultSortedFactor<>(m);
    assertFalse(s.isNormalized());
    assertEquals(0.4, s.p(0), 1e-10);
    assertEquals(0.6, s.p(1), 1e-10);
    assertEquals(1.0, s.p(2), 1e-10);
    assertEquals(hyp, s.hypotheses());
    assertEquals(arith, s.arithmetic());
    assertEquals(3, s.size());
    // and now the new stuff for a sorted factor ...
    assertEquals(1.0, s.value(0), 1e-10);
    assertEquals(0.6, s.value(1), 1e-10);
    assertEquals(0.4, s.value(2), 1e-10);
    assertEquals(2, s.hypothesis(0));
    assertEquals(1, s.hypothesis(1));
    assertEquals(0, s.hypothesis(2));
  }
}
