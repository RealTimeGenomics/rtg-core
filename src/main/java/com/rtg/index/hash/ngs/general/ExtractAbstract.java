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
package com.rtg.index.hash.ngs.general;

import java.io.IOException;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Deals with extracting masks.
 * Indel tweaking is done dynamically in the <code>maskIndel()</code> method.
 */
@TestClass({"com.rtg.index.hash.ngs.general.ExtractAbstractTest", "com.rtg.index.hash.ngs.general.SingleMaskTest"})
public abstract class ExtractAbstract extends IntegralAbstract {

  protected final SingleMask mSingleMask;

  ExtractAbstract(final SingleMask singleMask) {
    mSingleMask = singleMask;
  }

  /**
   * Called for each mask extraction.
   * Needed to accomodate the multiple masks that can result from allowance for indels.
   * @param bits final extracted bits using mask.
   * @throws IOException If an I/O error occurs
   */
  protected abstract void masked(final long bits) throws IOException;

  /**
   * Extract the bits from a single word.
   * @param v word bits are to be extracted from.
   * @return the extracted bits (right shifted aligned to position 0).
   */
  long mask(final long v) {
    assert mSingleMask.isFrozen();
    long res = 0;
    for (int i = 0; i < mSingleMask.size(); ++i) {
      final Skel sk = mSingleMask.subSkeleton(i);
      res = res | sk.mask(v);
    }
    return res;
  }

  /**
   * Extract the bits from a pair of words and pack the result into a single word.
   * @param v0 first word bits are to be extracted from.
   * @param v1 second word bits are to be extracted from.
   * @throws IOException If an I/O error occurs
   */
  void mask(final long v0, final long v1) throws IOException {
    assert mSingleMask.isFrozen();
    //System.err.println(mSingleMask.mWindowLength);
    final long res = (mask(v0) << mSingleMask.windowLength()) | mask(v1);
    masked(res);
  }

  //Used in recursive calls to ma
  private long mV0;

  private long mV1;

  /**
   * Extract the bits from a pair of words and pack the result into a single word.
   * Includes indel tweaking.
   * @param v0 first word bits are to be extracted from.
   * @param v1 second word bits are to be extracted from.
   * @throws IOException If an I/O error occurs
   */
  //TODO this code could probably do with some micro-optimization
  void maskIndel(final long v0, final long v1) throws IOException {
    assert mSingleMask.isFrozen();
    mV0 = v0;
    mV1 = v1;
//    System.err.println("v0:v1 " + Utils.toBits(v0, 12) + ":" + Utils.toBits(v1, 12) + " " + Utils.toCodes(v0, v1, 12));
    ma(0, 0, 0, mSingleMask.indels(), 0);
  }

  /**
   * @param r0 accumulated bits from <code>mV0</code> so far
   * @param r1 accumulated bits from <code>mV1</code> so far
   * @param i index of sub-mask
   * @param indels number of indels remaining to be used
   * @param tweak amount of tweak
   * @throws IOException If an I/O error occurs
   */
  private void ma(final long r0, final long r1, final int i, final int indels, final int tweak) throws IOException {
    assert indels >= 0;
    assert tweak <= mSingleMask.indels() * mSingleMask.indelLength() && -mSingleMask.indels()  * mSingleMask.indelLength() <= tweak;
    assert 0 <= i && i < mSingleMask.size();

    final Skel sk = mSingleMask.subSkeleton(i);
    final long s0 = r0 | sk.mask(mV0, tweak);
    final long s1 = r1 | sk.mask(mV1, tweak);
    final int j = i + 1;
    if (j == mSingleMask.size()) {
      //System.err.println("done i=" + i);
      //combine the two sets of extracted bits into a single word
      final long res = (s0 << mSingleMask.windowLength()) | s1;
      masked(res);
      return;
    }
    final int indelLength = mSingleMask.indelLength();
    final int indGap = indels * indelLength;
    final int min = Math.min(indGap, mSingleMask.gaps(i + 1));
    //System.err.println("indels=" + indels + " indelLength=" + indelLength + " indGap=" + indGap + " min=" + min);
    for (int t = -min; t < 0; ++t) {
      final int s = (-t + indelLength - 1) / indelLength;
      //System.err.println("neg i=" + i + " t=" + t + " s=" + s);
      ma(s0, s1, j, indels - s, tweak + t);
    }
    for (int t = 0; t <= indGap; ++t) {
      final int s = (t + indelLength - 1) / indelLength;
      //System.err.println("pos i=" + i + " t=" + t + " s=" + s);
      ma(s0, s1, j, indels - s, tweak + t);
    }
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("Extract");
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mSingleMask != null);
    return true;
  }
}
