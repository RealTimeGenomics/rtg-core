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

package com.rtg.variant.bayes.multisample.family;

import com.rtg.util.MathUtils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.bayes.Code;

/**
 * Class for producing a segregation score for use in VCF output.
 */
public class SegregationScore extends IntegralAbstract implements ISegregationScore {

  private final double[][] mLookup;
  private final int[][] mCounts;
  private final Code mCode;
  private final UniqueId mUid;
  private int mTotal = 0;

  /**
   * Constructor
   * @param code the code to use (should have range size of the number of alternative alleles plus one for the reference)
   * @param father the father hypothesis code
   * @param mother the mother hypothesis code
   */
  public SegregationScore(Code code, int father, int mother) {
    mCode = code;
    mUid = new UniqueId(mCode.rangeSize());
    final int fa0 = mUid.addId(mCode.a(father));
    assert fa0 == 0;
    final int fa1 = mUid.addId(mCode.bc(father));
    final int mo0 = mUid.addId(mCode.a(mother));
    final int mo1 = mUid.addId(mCode.bc(mother));

    mLookup = MendelianAlleleProbabilityDiploid.LOOKUP[fa1][mo0][mo1];
    mCounts = new int[mLookup.length][];
    for (int i = 0; i < mCounts.length; ++i) {
      mCounts[i] = new int[mLookup[i].length];
    }
  }

  @Override
  public void increment(int child, boolean diploid) {
    assert diploid;
    final int ch0 = mUid.id(mCode.a(child));
    if (ch0 < 0) {
      return;
    }
    final int ch1 = mUid.id(mCode.bc(child));
    if (ch1 < 0) {
      return;
    }
    mCounts[ch0][ch1]++;
    ++mTotal;
  }

  @Override
  public double lnProbability() {
    double lnp = 0.0;
    lnp += MathUtils.logFactorial(mTotal);
    for (int i = 0; i < mCounts.length; ++i) {
      for (int j = 0; j < mCounts[i].length; ++j) {
        final int count = mCounts[i][j];
        if (count > 0) {
          lnp -= MathUtils.logFactorial(count);
          lnp += mLookup[i][j] * count;
        }
      }
    }
    assert !Double.isNaN(lnp);
    return lnp;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mTotal >= 0);
    int t = 0;
    for (int[] mCount : mCounts) {
      for (final int count : mCount) {
        Exam.assertTrue(count >= 0);
        t += count;
      }
    }
    Exam.assertEquals(t, mTotal);
    return true;
  }

}
