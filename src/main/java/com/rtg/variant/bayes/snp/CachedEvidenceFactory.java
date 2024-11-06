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
