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

package com.rtg.util;

/**
 * Works out maximum aligner band width based on read length and a scaling factor.
 */
public class MaxShiftFactor {
  private final double mFactor;

  /**
   * Set scaling factor.
   * @param factor scaling factor.
   */
  public MaxShiftFactor(double factor) {
    if (factor < 0.0 || factor > 1) {
      throw new IllegalArgumentException("factor must be between 0 and 1: " + factor);
    }
    mFactor = factor;
  }

  /**
   * Returns the maximum shift for the give alignment threshold.
   * @param readLength the read length, to use for max shift threshold
   * @return maximum shift that is needed
   */
  public int calculateMaxShift(int readLength) {
    return Math.max((int) (readLength * mFactor), 7);
  }

  /**
   * Return the shift factor.
   * @return max shift factor.
   */
  public double getFactor() {
    return mFactor;
  }

  @Override
  public boolean equals(Object o) {
    return o instanceof MaxShiftFactor && Double.doubleToLongBits(((MaxShiftFactor) o).mFactor) == Double.doubleToLongBits(mFactor);
  }

  @Override
  public int hashCode() {
    return Double.valueOf(mFactor).hashCode() - 7;
  }

}
