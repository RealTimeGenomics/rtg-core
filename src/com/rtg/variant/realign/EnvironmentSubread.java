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
