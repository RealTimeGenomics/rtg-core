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

package com.rtg.variant.bayes.complex;

import com.rtg.sam.SamUtils;

/**
 */
public class MaxShiftCigarParser {

  private static final int MIN_MAX_SHIFT = 7; //calculated by Utils.calculateDefaultMaxShift(35), just is a good minimum width.

  private int mStartPos;
  private int mMaxShift;
  private int mSoftClipStartOffset;

  void parse(String cigar, int startPos) {
    int maxPosShift = 0;
    int maxNegShift = 0;
    mStartPos = startPos;
    mSoftClipStartOffset = 0;
    int count = 0;
    int currentOffset = 0;
    boolean isFirstAction = true;
    for (int i = 0; i < cigar.length(); ++i) {
      final char ch = cigar.charAt(i);
      if (Character.isDigit(ch)) {
        count = count * 10 + (ch - '0');
      } else {
        if (ch == SamUtils.CIGAR_INSERTION_INTO_REF) {
          currentOffset += count;
        } else if (ch == SamUtils.CIGAR_DELETION_FROM_REF) {
          currentOffset -= count;
        } else if (isFirstAction && ch == SamUtils.CIGAR_SOFT_CLIP) {
          mSoftClipStartOffset = count;
        }
        if (currentOffset > maxPosShift) {
          maxPosShift = currentOffset;
        } else if (currentOffset < maxNegShift) {
          maxNegShift = currentOffset;
        }
        count = 0;
        isFirstAction = false;
      }
    }

    mMaxShift = Math.max(MIN_MAX_SHIFT, (maxPosShift - maxNegShift) / 2 + 3);  //slop factor of 3 indels either side
    mStartPos = mStartPos - ((maxPosShift / 2) + (maxNegShift / 2));
//    System.err.println(cigar + "\t+" + maxPosShift + "\t" + maxNegShift + "\tms=" + mMaxShift + "\tsp=" + mStartPos);
  }

  int getMaxShift() {
    return mMaxShift;
  }

  int getStartPos() {
    return mStartPos;
  }

  int getSoftClipStartOffset() {
    return mSoftClipStartOffset;
  }
}
