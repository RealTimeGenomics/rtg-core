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

package com.rtg.reader;

/**
 * Computes a position to trim a read at via an algorithm very similar to:
 *
 * <code>argmax_x{\sum_{i=x+1}^l(INT-q_i)}</code> if <code>q_l&lt;INT</code> where l is the original read length.
 *
 */
public class BestSumReadTrimmer implements ReadTrimmer {

  private final int mQualityThreshold;

  /**
   * Construct a best sum read trimmer
   * @param qualityThreshold the threshold the sum the remainder of the read must be higher than
   */
  public BestSumReadTrimmer(int qualityThreshold) {
    mQualityThreshold = qualityThreshold;
  }

  @Override
  public int getTrimPosition(byte[] qualities, int length) {
    if (qualities.length == 0) {
      return 0;
    }
    int bestPos = length > qualities.length ? qualities.length : length;
    if (length > 0 && qualities[length - 1] < mQualityThreshold) {
      int bestSum = 0;
      int sum = 0;
      for (int i = length - 1; i >= 0; i--) {
        sum += mQualityThreshold - qualities[i];
        if (sum >= bestSum) {
          bestSum = sum;
          bestPos = i;
        }
      }
    }
    return bestPos;
  }
}
