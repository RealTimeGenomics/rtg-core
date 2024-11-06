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
 * Provides a probability distribution over a (haploid) set of hypotheses
 * that provide evidence for Bayesian estimation.
 * They will sum to approximately 1.0 and are in standard double format.
 */
public interface EvidenceInterface {

  /**
   * @param index selects the haploid hypothesis.
   * @return the probability that the specified hypothesis (index) is correct. (Ordinary double not in possibility arithmetic).
   */
  double probability(final int index);

  /**
   * @return the default probability that the specified evidence is correct in the absence of any information (including what the reference is). (Ordinary double not in possibility arithmetic).
   */
  double pe();

  /**
   * @return total probability that this is not equal to the read hypothesis. (Ordinary double not in possibility arithmetic).
   */
  double error();

  /**
   * @return probability that this read is not mapped to the correct position.  (Ordinary double not in possibility arithmetic).
   */
  double mapError();

  /**
   * @return the hypothesis corresponding to the sequence (often a single nucleotide)
   * being used to update the model, or <code>NOT_A_HYPOTHESIS</code> if the evidence
   * does not correspond to a single hypothesis.
   */
  int read();

  /**
   * number of bases to left of {@code read()} on read (i.e. to the left of the match)
   * @return number of bases to left of {@code read()} on read (i.e. to the left of the match)
   */
  int getReadBasesLeft();

  /**
   * number of bases to right of {@code read()} on read (i.e. to the right of the match)
   * @return number of bases to right of {@code read()} on read (i.e. to the right of the match)
   */
  int getReadBasesRight();

  /**
   * sets the value
   * @param readBasesLeft number of bases to left of {@code read()} on read (i.e. to the left of the match)
   */
  void setReadBasesLeft(int readBasesLeft);

  /**
   * sets the value
   * @param readBaseRight number of bases to right of {@code read()} on read (i.e. to the right of the match)
   */
  void setReadBasesRight(int readBaseRight);

  /**
   * @return if the evidence is mapped in forward frame
   */
  boolean isForward();

  /**
   * @return if the evidence is from a paired read
   */
  boolean isReadPaired();

  /**
   * @return if the evidence is from the first read of a pair
   */
  boolean isFirst();

  /**
   * @return if the evidence is from a properly mated paired read
   */
  boolean isMated();

  /**
   * @return if the evidence is from an unmapped arm
   */
  boolean isUnmapped();
}
