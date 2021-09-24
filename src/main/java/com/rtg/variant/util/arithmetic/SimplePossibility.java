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
