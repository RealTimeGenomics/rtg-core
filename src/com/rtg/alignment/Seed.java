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

package com.rtg.alignment;

import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Fast way of shifting in DNA to a shift register for use in EditDistance calculations.
 * Needs to deal with Ns hence the need for cunningness (a register with an N is invalid).
 */
public class Seed extends IntegralAbstract {

  private final int mSize;
  private final int mMask;

  /**
   * @param seedSize number of nucleotides in a seed.
   */
  Seed(final int seedSize) {
    assert seedSize >= 2;
    mSize = seedSize;
    final int seedBits = 2 * seedSize;
    final int mask0 = (1 << seedBits) - 1;
    mMask = mask0 | (mask0 << (Integer.SIZE - seedBits));
    integrity();
  }

  /**
   * @return the number of nucleotides in the seed.
   */
  public int size() {
    return mSize;
  }

  /**
   * @return initial value for a seed (everything is an N).
   */
  public int init() {
    return mMask;
  }

  /**
   * @param seed initial value of seed (all values are Ns).
   * @param nt nucleotide appended to seed.
   * @return new seed after nt appended and ofter oldest nt is deleted from start.
   */
  public int next(final int seed, final byte nt) {
    assert 0 <= nt && nt <= 4;
    return ((seed << 2) | (nt - 1)) & mMask;
  }

  /**
   * @param seed to be tested.
   * @return true iff the seed is valid, that is, contains no Ns.
   */
  public boolean isValid(final int seed) {
    return seed >= 0;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append("size=").append(mSize).append(" mask=").append(Utils.toBitsSep(mMask));
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mMask < 0);
    Exam.assertTrue(mSize * 4 + 2 <= Integer.SIZE);
    return true;
  }
}
