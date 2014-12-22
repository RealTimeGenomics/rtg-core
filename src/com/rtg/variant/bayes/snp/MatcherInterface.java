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

import com.rtg.variant.AbstractMachineErrorParams;
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
   * @param me error rates for a read group.
   * @param stateIndex the index number for the flag states
   */
  void match(int refPosition, int readBasesLeft, int readBasesRight, int readNt, int mapQ, int phred, AbstractMachineErrorParams me, int stateIndex);

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
   * @param isMated true if mapping from a paired read which is mated, false otherwise
   * @return the state index to pass to the match method
   */
  int getStateIndex(boolean isForward, boolean isReadPaired, boolean isMated);
}
