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
