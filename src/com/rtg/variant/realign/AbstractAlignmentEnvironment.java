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
