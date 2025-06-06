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
package com.rtg.index.hash;

import com.rtg.launcher.BuildParams;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class ExactHashFunction extends IntegralAbstract implements HashFunction {

  /** Number of bits in a long */
  private static final int LONG_BITS = 64;

  /** Number of codes to be included in the hash. */
  protected final int mWindowSize;

  /** A mask which includes just the bits needed to accomodate the window size. */
  protected final long mWindowMask;

  /** The number of bits needed to represent a single code. */
  protected final int mBits;

  private final long mCodeMask;

  /** The total number of bits needed to represent a complete window. */
  protected final int mWindowBits;

  private int mSoFar = 0;

  private final boolean mDualMode;

  protected long mHash;
  protected long mHashR;

  /**
   * @param buildParams holds the needed parameters.
   */
  public ExactHashFunction(final BuildParams buildParams) {
    this(buildParams.windowSize(), buildParams.sequences().mode().codeType().bits());
  }

  /**
   * @param buildParams holds the needed parameters.
   * @param dualMode whether to process forward and reverse simultaneously or not
   */
  public ExactHashFunction(final BuildParams buildParams, final boolean dualMode) {
    this(buildParams.windowSize(), buildParams.sequences().mode().codeType().bits(), dualMode);
  }

  /**
   * @param windowSize number of codes to be included in the window.
   * @param bits number of bits needed to represent a single code.
   */
  public ExactHashFunction(final int windowSize, final int bits) {
    this(windowSize, bits, false);
  }

  /**
   * @param windowSize number of codes to be included in the window.
   * @param bits number of bits needed to represent a single code.
   * @param dualMode whether to process forward and reverse simultaneously or not
   */
  public ExactHashFunction(final int windowSize, final int bits, final boolean dualMode) {
    mWindowSize = windowSize;
    mBits = bits;
    mWindowBits = mBits * mWindowSize;
    mWindowMask = (mWindowBits == LONG_BITS ? 0 : (1L << mWindowBits)) - 1;
    mCodeMask = (1L << mBits) - 1;
    mDualMode = dualMode;
    reset();
    integrity();
  }

  @Override
  public long hashStep(final byte code) {
    assert mWindowBits == LONG_BITS || (code >= 0 && code < (1 << mBits)) : "window bits=" + mWindowBits + " bits=" + mBits + " code=" + code;
    ++mSoFar;
    mHash = ((mHash << mBits) + code) & mWindowMask;
    assert mWindowBits == 64 || mHash >> mWindowBits == 0 : "hash=" + mHash + " windowBits=" + mWindowBits;
    if (mDualMode) {
      reverseHashStep(code);
    }
    return mHash;
  }

  long reverseHashStep(final byte code) {
    final int rc = 3 - code;
    mHashR = (mHashR >>> mBits) + ((long) rc << (mWindowBits - mBits));
    return mHashR;
  }

  /**
   * @return the hash
   */
  public final long hash() {
    return mHash;
  }

  /**
   * @return the reverse hash
   */
  public final long hashReverse() {
    return mHashR;
  }

  /**
   * Convert a hash to a String.
   * @param hash to be converted.
   * @param codes the character codes for the current encoding.
   * @return the string.
   */
  public String hashToSeq(final long hash, final char[] codes) {
    final char[] chars = new char[mWindowSize];
    long h = hash;
    for (int i = 0; i < mWindowSize; ++i) {
        final int x = (int) (h & mCodeMask);
        chars[mWindowSize - i - 1] = codes[x + 1];
        h = h >> mBits;
    }
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < mWindowSize; ++i) {
      sb.append(chars[i]);
    }
    return sb.toString();
  }

  /**
   * Check if the hash has been given sufficient steps to include a full window.
   * @return true iff the hash has been given sufficient steps to include a full window.
   */
  public final boolean isValid() {
    return mSoFar >= mWindowSize;
  }

  @Override
  public void reset() {
    mSoFar = 0;
    mHash = 0L;
    mHashR = 0L;
  }

  @Override
  public int getWindowSize() {
    return mWindowSize;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("SimpleHashFunction window size=").append(mWindowSize).append(" bits=").append(mBits).append(" mask=").append(mWindowMask).append(" hash=").append(mHash);
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mWindowSize >= 1 && mWindowSize <= LONG_BITS);
    Exam.assertTrue(mBits >= 1 && mBits <= LONG_BITS);
    Exam.assertTrue(mWindowBits >= 1 && mWindowBits <= LONG_BITS);
    Exam.assertEquals(0L, (mWindowMask + 1L) & mWindowMask);
    Exam.assertEquals(0L, (mWindowBits == LONG_BITS ? 0L : (1L << mWindowBits)) & mWindowMask);
    if ((1L << (mWindowBits - 1)) != ((1L << (mWindowBits - 1)) & mWindowMask)) {
      Exam.assertTrue(-1L + " " + this.toString(), false);
    }
    Exam.assertTrue(mSoFar >= 0);
    final int bitssofar = mSoFar * mBits;
    if (bitssofar < 64) {
      Exam.assertTrue(mHash  >>> bitssofar == 0);
    }
    return true;
  }
}

