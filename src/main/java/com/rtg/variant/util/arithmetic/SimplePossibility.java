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


/**
 * Implement as simple doubles.
 * This runs the risk of underflow which should be checked for and then calculations
 * done more slowly but carefully.
 * This class has no state only implements the methods
 * which manipulate doubles.
 */
public final class SimplePossibility implements PossibilityArithmetic {

  /** The only instance of this class that will ever be needed. */
  public static final PossibilityArithmetic SINGLETON = new SimplePossibility();

  private SimplePossibility() { }

  @Override
  public boolean isValidPoss(final double possibility) {
    return 0.0 <= possibility && possibility <= 1.0;
  }

  @Override
  public double prob2Poss(double value) {
    return value;
  }

  @Override
  public double poss2Prob(double possibility) {
    return possibility;
  }

  @Override
  public double ln2Poss(double value) {
    return Math.exp(value);
  }

  @Override
  public double poss2Ln(double possibility) {
    return Math.log(possibility);
  }
  @Override
  public double poss2Poss(double possibility, PossibilityArithmetic arithmetic) {
    return arithmetic.poss2Prob(possibility);
  }

  @Override
  public double zero() {
    return 0.0;
  }

  @Override
  public boolean isZero(double v) {
    return v == 0.0;
  }

  @Override
  public double one() {
    return 1.0;
  }

  @Override
  public double add(double v1, double v2) {
    return v1 + v2;
  }

  @Override
  public double multiply(double v1, double v2) {
    return v1 * v2;
  }

  @Override
  public double divide(double v1, double v2) {
    return v1 / v2;
  }

  @Override
  public double pow(double base, double exponent) {
    return Math.pow(base, exponent);
  }

  @Override
  public boolean underflow(double v) {
    return v < Double.MIN_NORMAL;
  }

  @Override
  public String toString() {
    return "SimplePossibility";
  }

  @Override
  public boolean gt(double v1, double v2) {
    return v1 > v2;
  }
}
