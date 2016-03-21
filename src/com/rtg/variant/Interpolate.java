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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.MathUtils;

/**
 * populates missing values in a list by drawing straight lines between adjacent values
 * Falls back on a default value (potentially per position) for leading/trailing missing values
 * @author kurt
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
    for (int i = minPos(); inBounds(i); i++) {
      if (!isMissing(i)) {
        for (int gap = prev + 1; gap < i; gap++) {
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
    for (int i = minPos(); inBounds(i); i++) {
      if (getValue(i) < 0) {
        setValue(i, values[i]);
      }
    }
  }
}
