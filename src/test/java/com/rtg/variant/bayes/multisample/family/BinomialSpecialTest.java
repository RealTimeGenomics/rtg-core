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

package com.rtg.variant.bayes.multisample.family;

import com.rtg.util.MathUtils;

import junit.framework.TestCase;

/**
 */
public class BinomialSpecialTest extends TestCase {

  public void test() {
    check(1, 1);
    check(2, 1);
    check(4, 2);
    check(4, 4);
    check(5, 5);
    check(6, 5);
    check(10, 5);
    check(100, 5);
    check(200, 5);
    check(BinomialSpecial.LENGTH - 2, 5);
    check(BinomialSpecial.LENGTH - 1, 5);
    check(BinomialSpecial.LENGTH, 5);
    check(BinomialSpecial.LENGTH + 1, 5);
  }

  private void check(int n, int a) {
    final double exp = MathUtils.logBinomial(n, a);
    final double actual = BinomialSpecial.logBinomial(n, a);
    //System.err.println(n + ":" + a + " " + Math.exp(exp) + ":" + Math.exp(actual));
    assertEquals(exp, actual, 1e-7);
  }
}
