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

/**
 * Provides a mapping from a hypothesis to an integer.
 */
public interface Code {

  /**
   * Extract the first part of the code from n.
   * @param n to be decoded.
   * @return first part of diploid code.
   */
  int a(int n);

  /**
   * Extract the second part of the code from n - only valid for diploid codes.
   * @param n to be decoded.
   * @return second part of diploid code.
   * @throws UnsupportedOperationException if diploid codes are not permitted.
   */
  int b(int n);

  /**
   * Extract the second part of the code from n - valid for all codes.
   * For a haploid code it returns the same value as <code>a(n)</code>.
   * @param n to be decoded.
   * @return second part of diploid code.
   */
  int bc(int n);

  /**
   * @param n code to be tested.
   * @return true iff n is a code for a homozygous hypothesis.
   */
  boolean homozygous(final int n);


  /**
   * @param hyp to be checked.
   * @return true iff <code>hyp</code> is in the valid range.
   */
  boolean valid(final int hyp);

  /**
   * @return the number of valid values in the code which ranges from 0..<code>size()</code> (exclusive).
   */
  int size();

  /**
   * @return the number of valid values returned by the <code>a()</code>, <code>b()</code> and <code>bc()</code> methods. They are all in the range
   * <code>[0..rangeSize())</code>.
   */
  int rangeSize();

  /**
   * Compute a code given a single homozygous hypothesis.
   * @param a the homozygous hypothesis.
   * @return a valid code.
   * @throws IllegalArgumentException if a does not specify a valid homozygous hypothesis.
   */
  int code(int a);

  /**
   * Compute a code given a pair of haploid hypothesis.
   * @param a a haploid hypothesis.
   * @param b a haploid hypothesis.
   * @return a valid code.
   * @throws IllegalArgumentException if a or b does not specify a valid haploid hypothesis.
   * @throws UnsupportedOperationException if the code does not support diploid codes and a &ne; b.
   */
  int code(int a, int b);
}
