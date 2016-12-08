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

import com.rtg.util.PortableRandom;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class InExactHashFunction extends IntegralAbstract implements
    HashFunction {

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
  public void reset() {
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

