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
