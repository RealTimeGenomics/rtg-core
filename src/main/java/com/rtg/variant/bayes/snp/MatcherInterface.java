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
 * Updates counts supplied from CIGAR format input.
 */
public interface MatcherInterface {

  /**
   * @param refPosition position on template (0 based)
   * @param readBasesLeft number of bases to left of <code>readNt</code> on read
   * @param readBasesRight number of bases to right of <code>readNt</code> on read
   * @param readNt read nucleotide in range 0=N, 1=A ... 5=D
   * @param mapQ phred score of probability that the read is not mapped to this position.
   * @param phred phred score of probability that the read nucleotide is not correct.
   * @param stateIndex the index number for the flag states
   */
  void match(int refPosition, int readBasesLeft, int readBasesRight, int readNt, int mapQ, int phred, int stateIndex);

  /**
   * Signal a piece of unmapped evidence.
   * @param refPosition position on template (0 based)
   */
  void unmapped(int refPosition);

  /**
   * @param refPosition position on template (0 based)
   * @param evid description of one piece of evidence.
   */
  void match(int refPosition, EvidenceInterface evid);

  /**
   * @param isForward true if mapped in forward frame, false otherwise
   * @param isReadPaired true if mapping is from a paired read, false otherwise
   * @param isFirst true if the mappings is single end of first of pair
   * @param isMated true if mapping from a paired read which is mated, false otherwise
   * @return the state index to pass to the match method
   */
  int getStateIndex(boolean isForward, boolean isReadPaired, boolean isFirst, boolean isMated);
}
