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

package com.rtg.variant.realign;

/**
 */
public interface AlignmentEnvironment {

  /**
   * Get the start position of the subsequence in the underlying sequence (usually a reference genome).
   * @return the start position of the read on the template (0 based).
   */
  int start();

  /**
   * Get the probability of the specified nucleotide
   * at the given index in the subsequence. These probabilities will typically
   * be computed from quality scores.
   * @param index position in subsequence (0 based - <code>start()</code> is 0).
   * @return the raw probability (not its log).
   */
  double quality(int index);

  /**
   * Get the base at specified position in the subsequence.
   * @param index position in sequence (0 based - <code>start()</code> is 0).
   * @return the nucleotide (0 = N).
   */
  byte base(int index);

  /**
   * Get a position which is guaranteed to be greater than the position of the last non-null nucleotide in the subsequence.
   * (0 based).
   * @return the length of the subsequence.
   */
  int subsequenceLength();

  /**
   * Return the length of the underlying template sequence.
   * @return template sequence length.
   */
  int templateLength();

  /**
   * Get whether the read is oriented forward. This is only necessary for CG reads where forward
   * corresponds to the overlap being on the left (lower template position).
   * @return true iff the the read is oriented forward.
   */
  boolean isInverted();

}
