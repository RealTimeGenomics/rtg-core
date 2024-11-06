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

