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
