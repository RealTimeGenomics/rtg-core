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
package com.rtg.complexity;

/**
 * Interface for classes that calculate sequence complexity scores.
 * Complexity scores are from SEG formula.
 *
 */
public interface SequenceComplexity {
  /**
   * Returns the minimum complexity from sliding window across sequence.
   *
   * @param sequence sequence to calculate complexity for
   * @return minimum complexity across sequence
   */
  double minComplexity(String sequence);

  /**
   * Returns the maximum complexity from sliding window across sequence.
   *
   * @param sequence sequence to calculate complexity for
   * @return maximum complexity across sequence
   */
  double maxComplexity(String sequence);

  /**
   * Returns the minimum complexity from sliding window across sequence.
   * The sequence is an encoded byte array as produced by <code>Utils.encodeArray</code>.
   *
   * @param sequence sequence to calculate complexity for
   * @return minimum complexity across sequence
   */
  double minComplexity(byte[] sequence);

  /**
   * Returns the maximum complexity from sliding window across sequence.
   * The sequence is an encoded byte array as produced by <code>Utils.encodeArray</code>.
   *
   * @param sequence sequence to calculate complexity for
   * @return maximum complexity across sequence
   */
  double maxComplexity(byte[] sequence);
}
