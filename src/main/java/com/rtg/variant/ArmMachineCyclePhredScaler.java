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

package com.rtg.variant;

import com.rtg.reader.Arm;

/**
 * Delegate phred scaling to separte quality scalers for each arm.
 */
public class ArmMachineCyclePhredScaler implements PhredScaler {
  PhredScaler mLeft;
  PhredScaler mRight;

  /**
   * @param left scaler for the left arm
   * @param right scaler for the right arm
   */
  public ArmMachineCyclePhredScaler(PhredScaler left, PhredScaler right) {
    mLeft = left;
    mRight = right;
  }

  @Override
  public int getScaledPhred(byte qual, int readPosition, Arm arm) {
    switch (arm) {
      case RIGHT:
        return mRight.getScaledPhred(qual, readPosition, arm);
      case LEFT:
      default:
        return mLeft.getScaledPhred(qual, readPosition, arm);
    }
  }
}
