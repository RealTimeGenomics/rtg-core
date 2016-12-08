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
import com.rtg.variant.util.VariantUtils;

/**
 */
public final class EvidenceQFactory implements CachedEvidenceFactory {

  private static final int NUM_STATES = 6;
  private static final int NUM_HYPOTHESIS = 4;
  private static final int MAX_MAPQ = 255;
  private static final int MAX_PHRED = 64;

  private static final EvidenceInterface[][][][] MEMO;
  static {
    MEMO = new EvidenceInterface[NUM_STATES][NUM_HYPOTHESIS][MAX_MAPQ][MAX_PHRED];
    for (int s = 0; s < NUM_STATES; ++s) {
      for (int i = 0; i < NUM_HYPOTHESIS; ++i) {
        for (int j = 0; j < MAX_MAPQ; ++j) {
          for (int k = 0; k < MAX_PHRED; ++k) {
            final boolean isForward = s >= 3;
            final boolean isReadPaired = (s % 3) >= 1;
            final boolean isMated = (s % 3) == 2;
            MEMO[s][i][j][k] = new EvidenceQ(DescriptionSnp.SINGLETON, i, VariantUtils.phredToProb(j), VariantUtils.phredToProb(k), isForward, isReadPaired, isMated, false);
          }
        }
      }
    }
  }

  /** Piece of evidence for an unmapped nucleotide. */
  private static final EvidenceInterface UNMAPPED = new EvidenceQ(DescriptionSnp.SINGLETON, 0, 1.0, 0, false, false, false, true); // most values irrelevant

  @Override
  public EvidenceInterface evidence(int readNt, int readBasesLeft, int readBaseRight, int mapQ, int phred, int stateIndex, int maxIndelLength, boolean isUnmapped) {
    if (isUnmapped) {
      return UNMAPPED;
    } else {
      final EvidenceInterface ret = MEMO[stateIndex][readNt][mapQ][phred >= MAX_PHRED ? MAX_PHRED - 1 : phred];
      ret.setReadBasesLeft(readBasesLeft);
      ret.setReadBasesRight(readBaseRight);
      return ret;
    }
  }

  /*
   * Possible states:
   * 0 - isForward=false, isReadPaired=false, isMated=false
   * 1 - isForward=false, isReadPaired=true, isMated=false
   * 2 - isForward=false, isReadPaired=true, isMated=true
   * 3 - isForward=true, isReadPaired=false, isMated=false
   * 4 - isForward=true, isReadPaired=true, isMated=false
   * 5 - isForward=true, isReadPaired=true, isMated=true
   * Illegal states:
   * isForward=true/false, isReadPaired=false, isMated = true
   */
  @Override
  public int getStateIndex(boolean isForward, boolean isReadPaired, boolean isMated) {
    int index = 0;
    if (isForward) {
      index += 3;
    }
    if (isReadPaired) {
      ++index;
      if (isMated) {
        ++index;
      }
    } else {
      assert !isMated;
    }
    return index;
  }

}
