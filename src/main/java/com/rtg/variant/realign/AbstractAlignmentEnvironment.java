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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.mode.DnaUtils;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
@TestClass(value = {"com.rtg.variant.realign.AlignmentEnvironmentReadTest"})
public abstract class AbstractAlignmentEnvironment extends IntegralAbstract implements AlignmentEnvironment {

  protected byte[] mRead;
  protected double[] mQuality;
  protected final int mStart;

  /**
   * @param start the position
   */
  public AbstractAlignmentEnvironment(int start) {
    mStart = start;
  }

  @Override
  public int start() {
    return mStart;
  }

  @Override
  public double quality(int index) {
    return mQuality[index];
  }

  @Override
  public byte base(int index) {
    return mRead[index];
  }

  @Override
  public int subsequenceLength() {
    return mRead.length;
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mRead);
    Exam.assertNotNull(mQuality);
    Exam.assertEquals(mRead.length, mQuality.length);
    return true;
  }

  @Override
  public String toString() {
    return "AlignmentEnvironment read=" + DnaUtils.bytesToSequenceIncCG(mRead) + " quality=" + Utils.realFormat(mQuality, 4) + " start=" + start();
  }

}
