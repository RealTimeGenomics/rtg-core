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
