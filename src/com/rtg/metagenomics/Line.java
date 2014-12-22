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
 */
public abstract class Line {

  /**
   * Compute the value at a specified position along the line.
   * @param delta position along line.
   * @return the value computed at that position.
   */
  public abstract double value(final double delta);

  /**
   * Get the value of the function and all derivatives that can be calculated.
   * @param delta position along line.
   * @return an array with the value at index 0 and any derivatives at order 1 ...
   */
  public double[] values(final double delta) {
    return new double[] {value(delta)};
  }

  /**
   * @return the highest order derivative that will be returned by <code>values</code> (&ge; 0).
   */
  public int derivativeOrder() {
    return 0;
  }
}
