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
 * A Double where the natural ordering is reversed (e.g. for sorting).
 */
public class ReverseDouble implements Comparable<ReverseDouble> {

  private final double mValue;

  /**
   * @param value value to reverse
   */
  public ReverseDouble(final double value) {
    mValue = value;
    if (Double.isNaN(value)) {
      throw new IllegalArgumentException();
    }
  }

  @Override
  public int compareTo(final ReverseDouble that) {
    return Double.compare(that.mValue, this.mValue);
  }

  /**
   * Returns the double value.
   * @return the double value.
   */
  public double doubleValue() {
    return mValue;
  }

  @Override
  public boolean equals(final Object obj) {
    return (obj instanceof ReverseDouble)
        && Double.doubleToLongBits(((ReverseDouble) obj).mValue) == Double.doubleToLongBits(this.mValue);
  }

  @Override
  public int hashCode() {
    final long lb = Double.doubleToRawLongBits(mValue);
    return (int) lb ^ (int) (lb >> 32);
  }

  @Override
  public String toString() {
    return Double.toString(mValue);
  }
}
