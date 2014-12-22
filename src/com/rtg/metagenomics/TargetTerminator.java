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

package com.rtg.metagenomics;

/**
 * Termination when solving for P-values. Stops when the estimate is signficantly smaller than the
 * distance to the parent L-value.
 */
class TargetTerminator implements Terminator {
  private final double mTerminationAbsolute;
  private final double mTerminationRatio;
  private final double mTarget;

  TargetTerminator(final double terminationAbsolute, final double terminationRatio, final double target) {
    mTerminationAbsolute = terminationAbsolute;
    mTerminationRatio = terminationRatio;
    mTarget = target;
  }

  @Override
  public boolean terminate(double estimate, double current) {
    assert current >= mTarget;
    return estimate < mTerminationRatio * (current - mTarget) || estimate < mTerminationAbsolute;
  }

}
