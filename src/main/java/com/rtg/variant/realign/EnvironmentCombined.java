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

import com.rtg.mode.DnaUtils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Combine two AlignmentEnvironments to form a combined environment ready for use by all paths.
 *
 */
public class EnvironmentCombined extends IntegralAbstract implements Environment {

  protected final AlignmentEnvironment mSamEnv;

  private final int mReadStart;

  protected final int mMaxShift;

  protected final AlignmentEnvironment mTemEnv;

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mSamEnv);
    return true;
  }

  /**
   * @param samEnv the SAM environment
   * @param zeroBasedStartPos the zero based read start pos
   * @param maxShift the maximum shift
   * @param templateEnv the template environment.
   */
  public EnvironmentCombined(AlignmentEnvironment samEnv, int zeroBasedStartPos, int maxShift, AlignmentEnvironment templateEnv) {
    mSamEnv = samEnv;
    mReadStart = zeroBasedStartPos;
    mMaxShift = maxShift;
    mTemEnv = templateEnv;
  }

  @Override
  public int absoluteTemplatePosition(int index) {
    return mReadStart + index;
  }

  @Override
  public int readLength() {
    return mSamEnv.subsequenceLength();
  }

  @Override
  public int maxShift() {
    return mMaxShift;
  }

  @Override
  public double quality(int index) {
    return mSamEnv.quality(index);
  }

  @Override
  public byte read(int index) {
    return mSamEnv.base(index);
  }

  @Override
  public byte template(int index) {
    return mTemEnv.base(index);
  }

  @Override
  public int templateLength() {
    return mTemEnv.templateLength();
  }

  @Override
  public void toString(StringBuilder sb) {
    super.toString(sb);
    sb.append("template after replace:").append(LS);
    final byte[] bs = new byte[templateLength()];
    int j = 0;
    for (int i = -mReadStart; i < templateLength() - mReadStart; ++i) {
      bs[j] = template(i);
      ++j;
    }
    sb.append(DnaUtils.bytesToSequenceIncCG(bs));
    sb.append(LS);
  }


}
