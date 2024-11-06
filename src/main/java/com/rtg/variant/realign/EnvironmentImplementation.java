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

import com.rtg.mode.DNA;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class EnvironmentImplementation extends IntegralAbstract implements Environment {

  protected static final double DEFAULT_QUALITY = 0.01; // for (CG) reads without quality.

  private static final char[] VALUE_CHARS = DNA.valueChars();

  protected final byte[] mTemplate;

  protected final byte[] mRead;

  /** The quality probabilities (generally 0 .. 0.25).  These are raw probabilities, not logs. */
  protected final double[] mQuality;

  private final int mStart;

  private final int mMaxShift;

  /**
   * @param maxShift maximum shift away from nominal start position (included so that <code>toString</code> is sensible).
   * @param template bases in template (coded as 0 = N ... 4 = T).
   * @param start nominal position the start of the read is mapped to.
   * @param read bases in read (coded as 0 = N ... 4 = T).
   * @param quality probability that the read base was incorrectly called by sequencing machine.
   */
  public EnvironmentImplementation(final int maxShift, final byte[] template, final int start, final byte[] read, final double[] quality) {
    super();
    mMaxShift = maxShift;
    mTemplate = template;
    //System.err.println("env start=" + start);
    mStart = start;
    mRead = read;
    mQuality = quality;
  }

  @Override
  public int readLength() {
    return mRead.length;
  }

  @Override
  public double quality(final int index) {
    if (mQuality == null) {
      return DEFAULT_QUALITY;
    }
    return mQuality[index];
  }

  @Override
  public byte read(final int index) {
    return mRead[index];
  }

  @Override
  public byte template(final int index) {
    final int i = mStart + index;
    if (i < 0 || i >= mTemplate.length) {
      return 0;
    }
    return mTemplate[i];
  }

  @Override
  public int absoluteTemplatePosition(final int index) {
    return mStart + index;
  }

  @Override
  public int maxShift() {
    return mMaxShift;
  }

  @Override
  public int templateLength() {
    return mTemplate.length;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("Environment").append(LS);
    final int end = readLength()  + maxShift();
    sb.append("Template [").append(absoluteTemplatePosition(0)).append("..").append(absoluteTemplatePosition(readLength())).append(")").append(LS);
    for (int i = -maxShift(); i < end; ++i) {
      final String indicator = i < 0 ? "-" : i >= readLength() ? "+" : " ";
      final int j = absoluteTemplatePosition(i);
      if (j >= 0 && j < templateLength()) {
        sb.append(indicator).append("[").append(i).append("]").append(VALUE_CHARS[template(i)]).append(LS);
      }
    }
    sb.append("Read").append(LS);
    for (int i = 0; i < readLength(); ++i) {
      sb.append("[").append(i).append("]").append(VALUE_CHARS[read(i)]).append(" ").append(Utils.realFormat(quality(i), 3)).append(LS);
    }
  }

  @Override
  public boolean integrity() {
    if (mQuality != null) {
      Exam.assertEquals(mRead.length, mQuality.length);
      Exam.assertTrue(readLength() <= mQuality.length);
    }
    Exam.assertTrue(readLength() <= mRead.length);
    Exam.assertTrue(templateLength() > 0);
    //Assert.assertTrue(0 <= templateStart() && templateStart() < templateLength());
    Exam.assertTrue(maxShift() >= 0);
    return true;
  }

}
