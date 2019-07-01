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
 * Factory for evidence objects with specified quality.
 */
public final class EvidenceQFactory implements CachedEvidenceFactory {

  private static final int NUM_STATES = 10;
  private static final int NUM_HYPOTHESIS = 4;
  private static final int MAX_MAPQ = 255;
  private static final int MAX_PHRED = 64;

  private static final EvidenceInterface[][][][] MEMO;
  static {
    MEMO = new EvidenceInterface[NUM_STATES][NUM_HYPOTHESIS][MAX_MAPQ][MAX_PHRED];
    for (int s = 0; s < NUM_STATES; ++s) {
      final boolean isForward = s >= 5;
      final int t = s % 5;
      final boolean isReadPaired = t > 0;
      final boolean isFirst = !isReadPaired || t > 2;
      final boolean isMated = isReadPaired && (t & 1) == 0;
      //System.out.println(s + " " + isForward + " " + isReadPaired + " " + isFirst + " " + isMated);
      for (int i = 0; i < NUM_HYPOTHESIS; ++i) {
        for (int j = 0; j < MAX_MAPQ; ++j) {
          for (int k = 0; k < MAX_PHRED; ++k) {
            MEMO[s][i][j][k] = new EvidenceQ(DescriptionSnp.SINGLETON, i, VariantUtils.phredToProb(j), VariantUtils.phredToProb(k), isForward, isReadPaired, isFirst, isMated, false);
          }
        }
      }
    }
  }

  /** Piece of evidence for an unmapped nucleotide. */
  private static final EvidenceInterface UNMAPPED = new EvidenceQ(DescriptionSnp.SINGLETON, 0, 1.0, 0, false, false, true, false, true); // most values irrelevant

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
   * 0 - isForward=false, isReadPaired=false, isFirst=true, isMated=false
   * 1 - isForward=false, isReadPaired=true, isFirst=false, isMated=false
   * 2 - isForward=false, isReadPaired=true, isFirst=false, isMated=true
   * 3 - isForward=false, isReadPaired=true, isFirst=true, isMated=false
   * 4 - isForward=false, isReadPaired=true, isFirst=true, isMated=true
   * 5 - isForward=true, isReadPaired=false, isFirst=true, isMated=false
   * 6 - isForward=true, isReadPaired=true, isFirst=false, isMated=false
   * 7 - isForward=true, isReadPaired=true, isFirst=false, isMated=true
   * 8 - isForward=true, isReadPaired=true, isFirst=true, isMated=false
   * 9 - isForward=true, isReadPaired=true, isFirst=true, isMated=true
   * Illegal states:
   * isForward=true/false, isReadPaired=false, isMated = true
   * isForward=true/false, isReadPaired=false, isFirst = false
   */
  @Override
  public int getStateIndex(boolean isForward, boolean isReadPaired, boolean isFirst, boolean isMated) {
    int index = 0;
    if (isForward) {
      index += 5;
    }
    if (isReadPaired) {
      ++index;
      if (isFirst) {
        index += 2;
      }
      if (isMated) {
        ++index;
      }
    } else {
      assert !isMated;
      assert isFirst;
    }
    return index;
  }

}
