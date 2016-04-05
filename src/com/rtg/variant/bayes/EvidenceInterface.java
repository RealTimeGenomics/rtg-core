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
   * @return if the evidence is from a properly mated paired read
   */
  boolean isMated();

  /**
   * @return if the evidence is from an unmapped arm
   */
  boolean isUnmapped();
}
