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
