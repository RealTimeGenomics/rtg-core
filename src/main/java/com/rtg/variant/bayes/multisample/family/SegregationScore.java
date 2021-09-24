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
