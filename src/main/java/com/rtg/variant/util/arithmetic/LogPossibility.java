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

package com.rtg.variant.util.arithmetic;

import com.rtg.variant.util.VariantUtils;

/**
 * Implement as natural logs.
 * Robust but slow. This class has no state only implements the methods
 * which manipulate doubles.
 */
public final class LogPossibility implements PossibilityArithmetic {

  /** The only instance of this class that will ever be needed. */
  public static final PossibilityArithmetic SINGLETON = new LogPossibility();

  private LogPossibility() { }

  @Override
  public boolean isValidPoss(final double possibility) {
    return possibility <= 0.0 && !Double.isNaN(possibility);
  }

  @Override
  public double prob2Poss(double value) {
    return Math.log(value);
  }

  @Override
  public double poss2Prob(double possibility) {
    return Math.exp(possibility);
  }

  @Override
  public double poss2Poss(double possibility, PossibilityArithmetic arithmetic) {
    return arithmetic.poss2Ln(possibility);
  }

  @Override
  public double ln2Poss(double value) {
    return value;
  }

  @Override
  public double poss2Ln(double possibility) {
    return possibility;
  }

  @Override
  public double zero() {
    return Double.NEGATIVE_INFINITY;
  }

  @Override
  public boolean isZero(double v) {
    return v == Double.NEGATIVE_INFINITY;
  }

  @Override
  public double one() {
    return 0.0;
  }

  @Override
  public double add(double v1, double v2) {
    return VariantUtils.logSum(v1, v2);
  }

  @Override
  public double multiply(double v1, double v2) {
    return v1 + v2;
  }

  @Override
  public double divide(double v1, double v2) {
    return v1 - v2;
  }

  @Override
  public double pow(double base, double exponent) {
    return exponent * base;
  }

  @Override
  public boolean underflow(double v) {
    return v == Double.NEGATIVE_INFINITY;
  }

  @Override
  public boolean gt(double v1, double v2) {
    return v1 > v2;
  }

  @Override
  public String toString() {
    return "LogPossibility";
  }
}
