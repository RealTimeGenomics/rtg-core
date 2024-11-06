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

/**
 */
public class BinomialAlleleBalanceTest extends AbstractAlleleBalanceTest {
  public void testAlleleBalanceZeroCoverage() {
    assertEquals(0.0, balance(new int[]{0, 0, 0, 0}, new double[]{0.0, 0.0, 0.0, 0.0}, findIndexByName("B", "B")));
  }
  public void testAlleleBalancePerfectBalance() {
    assertEquals(-3.370, balance(new int[]{9, 9, 0, 0}, new double[]{0.0, 0.0, 0.0, 0.0}, findIndexByName("B", "D")), 0.001);
  }

  public void testAlleleBalanceUnbalanced() {
    assertEquals(-12.477, balance(new int[]{9, 0, 0, 0}, new double[]{0.0, 0.0, 0.0, 0.0}, findIndexByName("B", "D")), 0.001);
  }

  public void testBinomialBalancedBetterThanUnbalanced() {
    final double balance = balance(new int[]{7, 7, 0, 0}, new double[]{0.1, 0.1, 0.0, 0.0}, findIndexByName("B", "D"));
    final double unbalanced = balance(new int[]{7, 4, 0, 0}, new double[]{0.1, 0.1, 0.0, 0.0}, findIndexByName("B", "D"));
    assertTrue(balance > unbalanced);
  }

  public void testBinomialGetsWorseAsBalanceDecreases() {
    final double unbalanced = balance(new int[]{7, 4, 0, 0}, new double[]{0.1, 0.1, 0.0, 0.0}, findIndexByName("B", "D"));
    final double worse = balance(new int[]{7, 3, 0, 0}, new double[]{0.1, 0.1, 0.0, 0.0}, findIndexByName("B", "D"));
    assertTrue(unbalanced > worse);
  }

  @Override
  protected BinomialAlleleBalance getAlleleBalanceProbability() {
    return new BinomialAlleleBalance(0.5);
  }
}
