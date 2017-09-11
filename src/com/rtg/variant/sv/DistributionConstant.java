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

package com.rtg.variant.sv;

import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;

/**
 * Constant.
 */
public class DistributionConstant extends Distribution {

  private final double mConstant;

  /**
   * @param lo low index of the distribution (inclusive).
   * @param hi high index of the distribution (exclusive).
   * @param constant value taken at all positions within the length.
   */
  public DistributionConstant(final int lo, final int hi, final double constant) {
    super(lo, hi);
    mConstant = constant;
    assert globalIntegrity();
  }

  @Override
  public Signal getSignalLn(SamCounts counts, String label) {
    return new SignalConstantLn(counts, this, label);
  }

  @Override
  protected double getValue(int index) {
    return mConstant;
  }

  /**
   * Returns the constant value for this distribution
   * @return the constant value
   */
  public double getConstant() {
    return mConstant;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append("Constant:").append(Utils.realFormat(mConstant, 4));
  }

  @Override
  public boolean integrity() {
    super.integrity();
    Exam.assertTrue(mConstant > 0.0 && Double.isFinite(mConstant));
    return true;
  }

}
