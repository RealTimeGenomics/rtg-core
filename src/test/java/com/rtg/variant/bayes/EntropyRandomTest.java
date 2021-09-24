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

import junit.framework.TestCase;

/**
 */
public class EntropyRandomTest extends TestCase {

  public void testEntropy() {
    checkEntropy(0.0, 0, 3);
    checkEntropy(Math.log(1.0 + 1.0 / 3.0), 1, 3);
    checkEntropy(2.0 * Math.log(2.0 + 1.0 / 5.0), 2, 5);
  }

  private void checkEntropy(final double exp, final int i, final int r) {
    assertEquals(exp, EntropyRandom.entropyLocal(i, r), 1e-10);
    assertEquals(exp, EntropyRandom.entropy(i, r), 1e-10);
  }

  public void testRandomPosterior() {
    checkRandomPosterior(0.0, 3);
    checkRandomPosterior(-2.0 * Math.log(2.0), 1, 1);
    //see spreadsheet
    checkRandomPosterior(-6.294491922, 3, 0, 2, 1);
  }

  private void checkRandomPosterior(final double exp, final int... counts) {
    int total = 0;
    for (final int c : counts) {
      total += c;
    }
    assertEquals(exp, EntropyRandom.randomPosteriorLocal(total, counts, 0.0), 1e-8);
    assertEquals(exp + 1.0, EntropyRandom.randomPosteriorLocal(total, counts, 1.0), 1e-8);
  }
}
