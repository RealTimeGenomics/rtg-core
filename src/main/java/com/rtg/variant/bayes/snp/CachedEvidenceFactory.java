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

package com.rtg.variant.bayes.snp;

import com.rtg.variant.bayes.EvidenceInterface;

/**
 * Factories for generating Probability objects depending on the underlying model
 * (does it support deletions or not).
 */
public interface CachedEvidenceFactory {

  /**
   * @param evidence the actual evidence, often a read nucleotide (0=A, 3=T, 4=D), but can be other things.
   * @param readBasesLeft number of bases to the left of <code>readNt</code> on read.
   * @param readBasesRight number of bases to the right of <code>readNt</code> on read.
   * @param mapQ the phred score for the probability that the read is not mapped to the current location.
   * @param phred the phred score for the probability that <code>readNt</code> is incorrectly called.
   * @param stateIndex the index number for the flag states
   * @param maxIndelLength  maximum length of an insert or delete as specified by cigar
   * @param isUnmapped true iff the evidence is unmapped
   * @return the probability.
   */
  EvidenceInterface evidence(int evidence, int readBasesLeft, int readBasesRight, int mapQ, int phred, int stateIndex, int maxIndelLength, boolean isUnmapped);

  /**
   * @param isForward true iff mapped in forward frame
   * @param isReadPaired true iff mapping is from a paired read
   * @param isFirst true iff the read in unpaired or the first of a pair
   * @param isMated true iff mapping from a paired read which is mated
   * @return the state index to pass to the evidence method
   */
  int getStateIndex(boolean isForward, boolean isReadPaired, boolean isFirst, boolean isMated);
}
