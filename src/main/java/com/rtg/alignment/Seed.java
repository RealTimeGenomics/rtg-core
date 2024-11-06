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
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Fast way of shifting in DNA to a shift register for use in edit distance calculations.
 * Needs to deal with Ns hence the need for cunningness (a register with an N is invalid).
 */
public class Seed extends IntegralAbstract {

  private final int mSize;
  private final int mMask; // Lower bits store the seed; upper bits store valid status (has an N been seen within the last mSize bases)

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
  public final boolean integrity() {
    Exam.assertTrue(mMask < 0);
    Exam.assertTrue(mSize * 4 + 2 <= Integer.SIZE);
    return true;
  }
}
