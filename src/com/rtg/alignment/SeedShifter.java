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
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Generates a succession of seeds from a sequence of DNA.
 */
public class SeedShifter extends IntegralAbstract {

  private final Seed mSeed;
  private final byte[] mBases;
  private final int mLength;
  private int mPosition;
  private int mValue;

  /**
   * @param seed used for calculating next seed.
   * @param bases array of nucleotides (assumed to be embedded in Ns before 0 and after length).
   * @param length number of nucleotides in <code>bases</code> that are in use
   * @param start position in DNA - may be negative or greater than length.
   */
  SeedShifter(Seed seed, byte[] bases, int length, int start) {
    mSeed = seed;
    mBases = bases;
    mLength = length;
    mPosition = start;
    mValue = mSeed.init();
  }

  /**
   * @return seed with next nucleotide sifted in.
   */
  public int next() {
    final byte nt;
    if (mPosition < 0 || mPosition >= mLength) {
      nt = 0;
    } else {
      nt = mBases[mPosition];
    }
    ++mPosition;
    mValue = mSeed.next(mValue, nt);
    return mValue;
  }

  /**
   * @return the position of the first nucleotide included in the shift.
   */
  public int position() {
    return mPosition - mSeed.size();
  }

  /**
   * @return true iff the seed is valid (contains no Ns).
   */
  public boolean isValid() {
    return mSeed.isValid(mValue);
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append("position=").append(mPosition).append(" value=").append(Utils.toBitsSep(mValue));
  }

  @Override
  public boolean integrity() {
    return true;
  }

}
