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

package com.rtg.variant.sv;

import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;

/**
 * Constant.
 */
public final class DistributionConstant extends Distribution {

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
