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

package com.rtg.variant;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.MathUtils;

/**
 * populates missing values in a list by drawing straight lines between adjacent values
 * Falls back on a default value (potentially per position) for leading/trailing missing values
 */
@TestClass({"com.rtg.variant.Interpolate2dArrayTest", "com.rtg.variant.InterpolateArrayTest"})
interface Interpolate {
  /**
   * @param pos index to fetch
   * @return the value at index {@code pos}
   */
  int getValue(int pos);

  /**
   * @param pos index to set
   * @param value the value to store
   */
  void setValue(int pos, int value);

  /**
   * @return the smallest index to interpolate
   */
  int minPos();

  /**
   * @param pos index to check
   * @return true if the specified position is within bounds
   */
  boolean inBounds(int pos);

  /**
   * @param pos index to check
   * @return true if the value at {@code pos} represents a missing value
   */
  boolean isMissing(int pos);

  /**
   *  Fill missing values in your structure by taking points on the straight line between present values
   */
  default void interpolate() {
    boolean foundFirstValue = false;
    int prev = minPos() - 1;
    for (int i = minPos(); inBounds(i); ++i) {
      if (!isMissing(i)) {
        for (int gap = prev + 1; gap < i; ++gap) {
          if (foundFirstValue) {
            final int gapSize = i - prev;
            final int gapValue = getValue(i) - getValue(prev);
            final int step = gap - prev;
            final int qual = getValue(prev) + (int) MathUtils.round((double) gapValue / gapSize * step);
            setValue(gap, qual);
          }
        }
        prev = i;
        foundFirstValue = true;
      }
    }
  }

  default void fill(int[] values) {
    for (int i = minPos(); inBounds(i); ++i) {
      if (getValue(i) < 0) {
        setValue(i, values[i]);
      }
    }
  }
}
