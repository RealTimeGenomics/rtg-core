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

package com.rtg.variant.realign;

import com.rtg.util.integrity.Exam;

/**
 * An environment that uses a sub-read of the read array.
 */
public class EnvironmentSubread extends EnvironmentImplementation {

  private final int mReadStart;

  private final int mReadEnd;

  /**
   * @param maxShift maximum shift away from nominal start position (included so that <code>toString</code> is sensible).
   * @param template bases in template (coded as 0 = N ... 4 = T).
   * @param start nominal position the start of the read is mapped to.
   * @param read bases in read (coded as 0 = N ... 4 = T).
   * @param quality probability that the read base was incorrectly called by sequencing machine.
   */
  public EnvironmentSubread(final int maxShift, final byte[] template, final int start, final byte[] read, final double[] quality) {
    this(maxShift, template, start, read, 0, read.length, quality);
  }

  /**
   * @param maxShift maximum shift away from nominal start position (included so that <code>toString</code> is sensible).
   * @param template bases in template (coded as 0 = N ... 4 = T).
   * @param start nominal position the start of the read is mapped to.
   * @param read bases in read (coded as 0 = N ... 4 = T).
   * @param readStart start position of read (0 based).
   * @param readEnd end position of read (0 based exclusive).
   * @param quality probability that the read base was incorrectly called by sequencing machine.
   */
  public EnvironmentSubread(final int maxShift, final byte[] template, final int start, final byte[] read, final int readStart, final int readEnd, final double[] quality) {
    super(maxShift, template, start, read, quality);
    mReadStart = readStart;
    mReadEnd = readEnd;
  }

  @Override
  public int readLength() {
    return mReadEnd - mReadStart;
  }

  @Override
  public double quality(final int index) {
    if (mQuality == null) {
      return DEFAULT_QUALITY;
    }
    final int ix = index + mReadStart;
    assert 0 <= ix && ix < mReadEnd;
    return mQuality[ix];
  }

  @Override
  public byte read(final int index) {
    final int ix = index + mReadStart;
    assert 0 <= ix && ix < mReadEnd;
    return mRead[ix];
  }

  @Override
  public boolean integrity() {
    super.integrity();
    Exam.assertTrue(0 <= mReadStart && mReadStart <= mReadEnd && mReadEnd <= super.readLength());
    return true;
  }

}
