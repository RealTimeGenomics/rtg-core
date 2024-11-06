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
package com.rtg.variant.bayes.multisample;

import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;

/**
 * This version looks for how much of the subsequence can be explained by a single
 * n-mer obtained from the center of the subsequence
 */
public final class SingleNMerRepeatMeasurer implements RepeatMeasurer {

  private static final int MAX_MER_LENGTH = GlobalFlags.getIntegerValue(CoreGlobalFlags.COMPLEX_REGION_SIMPLE_REPEAT_LIMIT);

  private final byte[] mReferenceNts;

  private final int mMaxMerLength;

  /**
   * Constructor
   * @param referenceNts reference sequence bytes
   */
  public SingleNMerRepeatMeasurer(byte[] referenceNts) {
    this(referenceNts, MAX_MER_LENGTH);
  }

  /**
   * Constructor
   * @param referenceNts reference sequence bytes
   * @param maxMerLength length of longest n-mer simple repeat type
   */
  public SingleNMerRepeatMeasurer(byte[] referenceNts, int maxMerLength) {
    mReferenceNts = referenceNts;
    mMaxMerLength = maxMerLength;
  }

  @Override
  public int measureRepeats(int positionA, int positionB) {
    return measureRepeats(positionA, positionB, 0);
  }

  @Override
  public int measureRepeats(int positionA, int positionB, int repeatHint) {
    final int startPos = Math.max(positionA, 0);
    final int endPos = Math.min(positionB, mReferenceNts.length);
    final int gapLength = endPos - startPos;
    final int midpoint = startPos + gapLength / 2;
    final int maxMerLength = repeatHint == 0 ? mMaxMerLength : repeatHint + 1;
    int repeatTotal = 0;
    for (int simpleLength = 1; simpleLength <= maxMerLength; ++simpleLength) {
      if (simpleLength > gapLength) { // Dont try to look for repeat units larger than the gap size
        break;
      }

      final int repeatStart = midpoint - simpleLength / 2;
      int intervalMatches = simpleLength; // Count of bases within the interval that are part of the repeat unit
      int extervalMatches = 0;            // Count of bases that continue the repeat unit outside the interval

      // Scan towards the end
      int repeatPos = repeatStart;
      int copyPos = repeatStart + simpleLength;
      while ((copyPos < endPos) && (mReferenceNts[repeatPos] != 0) && (mReferenceNts[repeatPos++] == mReferenceNts[copyPos++])) {
        ++intervalMatches;
      }
      if (copyPos == endPos) {
        // See if the repeat continues right past the end pos
        while ((repeatPos < endPos) && (copyPos < mReferenceNts.length) && (mReferenceNts[repeatPos++] == mReferenceNts[copyPos++])) {
          ++extervalMatches;
        }
      }

      // Scan towards the start
      repeatPos = repeatStart + simpleLength - 1;
      copyPos = repeatStart - 1;
      while ((copyPos >= startPos) && (mReferenceNts[repeatPos] != 0) && (mReferenceNts[repeatPos--] == mReferenceNts[copyPos--])) {
        ++intervalMatches;
      }
      if (copyPos < startPos) {
        // See if the repeat continues left of start pos
        while ((repeatPos >= endPos) && (copyPos >= 0) && (mReferenceNts[repeatPos--] == mReferenceNts[copyPos--])) {
          ++extervalMatches;
        }
      }

      final int targetBases = simpleLength + simpleLength / 2; // Only recognize a pattern if there is a repeat of at least half the n-mer
      if (intervalMatches + extervalMatches > targetBases && intervalMatches > repeatTotal) {
        repeatTotal = intervalMatches;
        //System.err.println("Total repeats: " + repeatTotal + " via r" + simpleLength + " " + DnaUtils.bytesToSequenceIncCG(mReferenceNts, repeatStart, simpleLength));
        if (repeatTotal == gapLength) { // Already explained the whole gap (e.g. for aaaaaaaaaaaa no need to try 2-mer/3-mer/4-mer, etc)
          break;
        }
      }
    }
    return repeatTotal;
  }
}
