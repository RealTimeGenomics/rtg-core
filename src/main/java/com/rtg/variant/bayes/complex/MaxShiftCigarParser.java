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
