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
