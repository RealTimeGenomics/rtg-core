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
package com.rtg.pairedend;


/**
 * Extends Hit Info class for extra data necessary for unfiltered paired end
 */
public class UnfilteredHitInfo extends AbstractHitInfo<UnfilteredHitInfo> {

  private boolean mMatedOk = false;
  private boolean mUnmatedOk = false;

  /**
   */
  @Override
  public void setValues(boolean first, boolean reverseComplement, int readId, int templateStart) {
    super.setValues(first, reverseComplement, readId, templateStart);
    mMatedOk = false;
    mUnmatedOk = false;
  }

  /**
   * @param value true if hit has passed the mated threshold.
   */
  public void setMatedOk(boolean value) {
    mMatedOk = value;
  }

  /**
   * @return true if this hit has passed the mated threshold.
   */
  public boolean getMatedOk() {
    return mMatedOk;
  }

  /**
   * @param value true if hit has passed the unmated threshold.
   */
  public void setUnmatedOk(boolean value) {
    mUnmatedOk = value;
  }

  /**
   * @return true if this hit has passed the unmated threshold.
   */
  public boolean getUnmatedOk() {
    return mUnmatedOk;
  }
}
