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
 * A possibility is like a probability where you don't have to know how it is represented
 * (other than it is packed into a double). This interface is for arithmetic operations on
 * these possibilities. The intention is to deal with numbers with very large ranges of
 * exponents.
 */
public interface PossibilityArithmetic {

  /**
   * @param possibility an internal Possibility value.
   * @return true iff it is valid.
   */
  boolean isValidPoss(final double possibility);


  /**
   * Convert from a potentially different arithmetic space into this one
   * @param possibility a possibility in the source space
   * @param arithmetic the source arithmetic space
   * @return a possibility in this space
   */
  double poss2Poss(double possibility, PossibilityArithmetic arithmetic);

  /**
   * @param value any valid probability, from 0 to 1.0
   * @return an internal Possibility value.
   */
  double prob2Poss(double value);

  /**
   * @param possibility an internal Possibility value.
   * @return a valid probability, from 0 to 1.0
   */
  double poss2Prob(double possibility);

  /**
   * @param value the natural log of a probability.
   * @return an internal Possibility value.
   */
  double ln2Poss(double value);

  /**
   * @param possibility an internal Possibility value.
   * @return the natural log of a probability.
   */
  double poss2Ln(double possibility);

  /**
   * @return a possibility representing the probability 0.0.
   */
  double zero();

  /**
   * @param v value to be tested.
   * @return true iff v is zero.
   */
  boolean isZero(double v);

  /**
   * @return a possibility representing the probability 1.0.
   */
  double one();

  /**
   * @param v1 an internal Possibility value.
   * @param v2 an internal Possibility value.
   * @return the sum of <code>v1</code> and <code>v2</code> as a possibility.
   */
  double add(double v1, double v2);

  /**
   * @param v1 an internal Possibility value.
   * @param v2 an internal Possibility value.
   * @return the product of <code>v1</code> and <code>v2</code> as a possibility.
   */
  double multiply(double v1, double v2);

  /**
   * @param v1 an internal Possibility value.
   * @param v2 an internal Possibility value.
   * @return <code>v1</code> divided by <code>v2</code> as a possibility.
   */
  double divide(double v1, double v2);

  /**
   * @param base arithmetic value
   * @param exponent exponent in natural form (i.e. not in arithmetic)
   * @return <code>base^exponent</code>
   */ 
  double pow(double base, double exponent);

  /**
   * Check if the result of a series of multiplications has underflowed.
   * @param v the value to be checked.
   * @return true iff v has underflowed.
   */
  boolean underflow(double v);

  /**
   * @param v1 first value to be checked.
   * @param v2 second value to be checked.
   * @return true iff <code>v1</code> is greater than <code>v2</code> in probability space.
   */
  boolean gt(double v1, double v2);
  
}
