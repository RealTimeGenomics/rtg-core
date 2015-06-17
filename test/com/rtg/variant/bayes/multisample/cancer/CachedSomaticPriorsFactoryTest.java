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

package com.rtg.variant.bayes.multisample.cancer;

import java.io.IOException;
import java.util.Arrays;

import com.rtg.util.InvalidParamsException;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.MockHypotheses;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class CachedSomaticPriorsFactoryTest extends TestCase {

  public void test() throws InvalidParamsException, IOException {
    final DescriptionCommon desc = new DescriptionCommon("", "A", "AA");
    final Hypotheses<DescriptionCommon> hyp = new MockHypotheses<>(desc, SimplePossibility.SINGLETON, true, new double[] {0.0, 0.0, 0.0}, 0);
    final CachedSomaticPriorsFactory<DescriptionCommon> cache = new CachedSomaticPriorsFactory<>(hyp, 0.0);
    final double[][] qfc = cache.somaticQ(0.5);
    final double[][] qf = new SomaticPriorsFactory<>(hyp, 0.0).somaticQ(0.5); // This value is not binned since it is a power of 2
    assertTrue(Arrays.deepEquals(qfc, qf));
    assertTrue(qfc == cache.somaticQ(0.5));
    assertTrue(qfc == cache.somaticQ(0.5 + Double.MIN_NORMAL));
    assertFalse(qfc == cache.somaticQ(0.25));
    // Check that extreme mu doesn't cause an exception
    cache.somaticQ(1);
    cache.somaticQ(Double.MIN_VALUE);
    cache.somaticQ(0.0);
  }

  public void testBounds() {
    assertEquals(-1023, Math.getExponent(0.0));
    assertEquals(-1023, Math.getExponent(Double.MIN_VALUE));
    assertEquals(-1022, Math.getExponent(Double.MIN_NORMAL));
    assertEquals(0, Math.getExponent(1.0));
  }

}
