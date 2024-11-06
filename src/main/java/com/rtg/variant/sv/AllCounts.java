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

package com.rtg.variant.sv;

import static com.rtg.util.StringUtils.TAB;

import java.io.PrintStream;

import com.reeltwo.jumble.annotations.JumbleIgnore;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
class AllCounts extends IntegralAbstract {

  private static final int PROPER_LEFT = 0;
  private static final int DISCORD_LEFT = 1;
  private static final int UNMATED_LEFT = 2;
  private static final int PROPER_RIGHT = 3;
  private static final int DISCORD_RIGHT = 4;
  private static final int UNMATED_RIGHT = 5;
  private static final int UNPAIRED = 6;
  private static final int LENGTH = 7;

  final SamArray[] mCounts = new SamArray[LENGTH];

  private final int mLength;

  AllCounts(int length) {
    mLength = length;
    for (int i = 0; i < LENGTH; ++i) {
      mCounts[i] = new SamArray(length);
    }
  }

  SamArray properLeft() {
    return mCounts[PROPER_LEFT];
  }

  SamArray discordantLeft() {
    return mCounts[DISCORD_LEFT];
  }

  SamArray unmatedLeft() {
    return mCounts[UNMATED_LEFT];
  }

  SamArray properRight() {
    return mCounts[PROPER_RIGHT];
  }

  SamArray discordantRight() {
    return mCounts[DISCORD_RIGHT];
  }

  SamArray unmatedRight() {
    return mCounts[UNMATED_RIGHT];
  }

  SamArray unpaired() {
    return mCounts[UNPAIRED];
  }

  void plot(PrintStream ps) {
    for (int i = 0; i < mLength; ++i) {
      ps.print(i);
      for (int j = 0; j < LENGTH; ++j) {
        ps.print(TAB + Utils.realFormat(mCounts[j].count(0, i), 4));
      }
      ps.println();
    }
  }

  @Override
  public boolean integrity() {
    Exam.assertEquals(7, mCounts.length);
    Exam.assertTrue(mLength >= 0);
    for (int i = 0; i < LENGTH; ++i) {
      Exam.assertEquals(mLength, mCounts[i].length());
    }
    return true;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append("AllCounts:").append(mLength);
  }

  /**
   * Reverse the counts around their center and swap
   * left and right arms. This gives the other breakpoint
   * for asymmetric signals.
   * @return the reversed signals.
   */
  AllCounts reverse() {
    final AllCounts rev = new AllCounts(mLength);
    for (int i = 0; i < mLength; ++i) {
      final int j = mLength - i - 1;
      rev.properLeft().increment(j, properRight().count(0, i));
      //System.err.println("i=" + i + " j=" + j + " count=" + properRight().count(0, i));
      rev.discordantLeft().increment(j, discordantRight().count(0, i));
      rev.unmatedLeft().increment(j, unmatedRight().count(0, i));
      rev.properRight().increment(j, properLeft().count(0, i));
      rev.discordantRight().increment(j, discordantLeft().count(0, i));
      rev.unmatedRight().increment(j, unmatedLeft().count(0, i));
      rev.unpaired().increment(j, unpaired().count(0, i));
    }
    return rev;
  }

  /**
   * Completes the right hand sides as the left hand sides reflected around a pivot.
   * Used for testing.
   * @param pivot point around which the arrays are swapped.
   * @param offset to apply after swapping.
   */
  @JumbleIgnore
  void completeRight(final int pivot, final int offset) {
    mCounts[PROPER_RIGHT] = mCounts[PROPER_LEFT].reverse(pivot, offset);
    mCounts[DISCORD_RIGHT] = mCounts[DISCORD_LEFT].reverse(pivot, offset);
    mCounts[UNMATED_RIGHT] = mCounts[UNMATED_LEFT].reverse(pivot, offset);
  }
}
