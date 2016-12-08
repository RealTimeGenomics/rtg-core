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
