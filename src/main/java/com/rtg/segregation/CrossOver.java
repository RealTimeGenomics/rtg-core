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

package com.rtg.segregation;


import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public final class CrossOver extends IntegralAbstract {

  private final boolean mFather;

  private final int mChild;

  private final int mCountBefore;

  private final int mCountAfter;

  private final int mLength;

  private final PatternArray mConsistent;

  /**
   * @param length number of children in family.
   * @param father true iff crossover is from father.
   * @param child the child the crossover occurred in.
   * @param countBefore of one arbitrary side of the phasing before the cross-over (be careful this is tricky to use).
   * @param countAfter of the same side of the phasing after the cross-over (be careful this is tricky to use).
   * @param consistent pattern after the crossover that has consistent labeling with the pattern before the crossover.
   */
  CrossOver(int length, boolean father, int child, int countBefore, int countAfter, PatternArray consistent) {
    mLength = length;
    mFather = father;
    mChild = child;
    mCountBefore = countBefore;
    mCountAfter = countAfter;
    mConsistent = consistent;
    integrity();
  }

  /**
   * Compute the next levels for father and mother and keep as compatible as possible given the previous levels (there is ambiguity because
   * we dont know the absolute phasing).
   * @param faPrev the previous level for the father from the search output logic.
   * @param moPrev the previous level for the mother from the search output logic.
   * @return new level counts (father then mother).
   */
  int[] searchLevel(final int faPrev, final int moPrev) {
    final int faAfter;
    final int moAfter;
    assert 0 <= faPrev && faPrev <= mLength;
    assert 0 <= moPrev && moPrev <= mLength;
    if (mFather) {
      if (mCountBefore == faPrev) {
        faAfter = mCountAfter;
      } else if (mCountBefore == mLength - faPrev) {
        faAfter = mLength - mCountAfter;
      } else {
        throw new RuntimeException(); //"before=" + mCountBefore + " after=" + mCountAfter + " faPrev=" + faPrev + " moPrev=" + moPrev);
      }
      moAfter = moPrev;
    } else {
      if (mCountBefore == moPrev) {
        moAfter = mCountAfter;
      } else if (mCountBefore == mLength - moPrev) {
        moAfter = mLength - mCountAfter;
      } else {
        throw new RuntimeException(); //"before=" + mCountBefore + " after=" + mCountAfter + " faPrev=" + faPrev + " moPrev=" + moPrev);
      }
      faAfter = faPrev;
    }
    return new int[] {faAfter, moAfter};
  }

  PatternArray pattern() {
    return mConsistent;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append(mFather ? "fa" : "mo");
    sb.append(" ").append(mChild);
    sb.append(" ").append(mCountBefore);
    sb.append(" ").append(mCountAfter);
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mChild >= 0);
    Exam.assertTrue(0 <= mCountBefore && mCountBefore <= mLength);
    Exam.assertTrue(0 <= mCountAfter && mCountAfter <= mLength);
    Exam.assertTrue(Math.abs(mCountAfter - mCountBefore) == 1);
    Exam.assertNotNull(mConsistent);
    return true;
  }
}
