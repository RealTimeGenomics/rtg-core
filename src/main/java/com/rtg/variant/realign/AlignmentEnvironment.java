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
