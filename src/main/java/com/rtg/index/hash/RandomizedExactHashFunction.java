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
package com.rtg.index.hash;

/**
 * The aim of this is to give an exact hash that has its high order bits well randomized
 * so that they can effectively be used as indices.
 * This is not perfect in particular the bits of the randomized version will depend only on
 * the bits to the right in the original. But it is fast and simple.
 */
public class RandomizedExactHashFunction extends ExactHashFunction {

  private final long mPrime;

  /**
   * @param windowSize number of codes to be included in the window.
   * @param bits number of bits needed to represent a single code.
   */
  public RandomizedExactHashFunction(final int windowSize, final int bits) {
    super(windowSize, bits);
    mPrime = PrimeUtils.prime(windowSize);
  }

  @Override
  public long hashStep(final byte code) {
    final long ha = super.hashStep(code);
    final long res = (ha * mPrime) & mWindowMask;
    assert mWindowBits == 64 || res >> mWindowBits == 0;
    //System.err.println("hashStep code=" + code + " ha=" + ha + " res=" + res);
    return res;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("RandomizedExactHashFunction window size=").append(mWindowSize).append(" bits=").append(mBits).append(" mask=").append(mWindowMask).append(" hash=").append(mHash).append(" prime=").append(mPrime);
  }

}

