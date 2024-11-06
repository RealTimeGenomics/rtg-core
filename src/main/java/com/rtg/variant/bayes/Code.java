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
