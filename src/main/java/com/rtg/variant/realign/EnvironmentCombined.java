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
