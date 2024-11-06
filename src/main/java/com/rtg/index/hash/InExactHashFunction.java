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

import com.rtg.util.PortableRandom;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class InExactHashFunction extends IntegralAbstract implements HashFunction {

  /** Randomly generated arrays used to compute <code>irvineHash</code> codes */
  private static final long[] HASH_BLOCKS;
  static {
    HASH_BLOCKS = new long[256];
    final PortableRandom r = new PortableRandom(1); // use same seed for deterministic behavior
    for (int i = 0; i < 256; ++i) {
      HASH_BLOCKS[i] = r.nextLong();
    }
  }

  /** Number of codes to be included in the hash. */
  private final int mWindowSize;

  private int mI;

  private long mHash;

  /**
   * @param windowSize number of codes to be included in the window.
   */
  public InExactHashFunction(final int windowSize) {
    mWindowSize = windowSize;
    reset();
  }

  @Override
  public long hashStep(final byte code) {
    assert mI < mWindowSize;
    mHash = Long.rotateLeft(mHash, 1) ^ HASH_BLOCKS[(code + mI++) & 0xFF];
    return mHash;
  }

  @Override
  public final void reset() {
    mHash = 0L;
    mI = 0;
  }

  @Override
  public int getWindowSize() {
    return mWindowSize;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("IrvineHashFunction window size=").append(mWindowSize).append(" hash=").append(mHash).append(" i=").append(mI);
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mWindowSize >= 1);
    Exam.assertTrue(mI >= 0 && mI <= mWindowSize);
    return true;
  }

}

