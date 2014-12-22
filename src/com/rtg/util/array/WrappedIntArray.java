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
package com.rtg.util.array;

import java.util.Arrays;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Wraps an integer array so that it is immutable.
 */
public class WrappedIntArray extends IntegralAbstract implements ImmutableIntArray {

  private final int[] mArray;

  /**
   * @param array to be wrapped.
   */
  public WrappedIntArray(final int[] array) {
    mArray = array.clone();
  }

  /**
   * @param array to be wrapped.
   */
  public WrappedIntArray(final long[] array) {
    mArray = new int[array.length];
    for (int i = 0; i < array.length; i++) {
      mArray[i] = (int) array[i];
    }
  }

  /**
   * Construct lengths array where lengths for left and right arms alternate.
   * @param left lengths  of reads on left arm.
   * @param right lengths of reads on right arm.
   */
  public WrappedIntArray(final int[] left, final int[]  right) {
    assert left.length == right.length;
    final int size = left.length + right.length;
    mArray = new int[size];
    for (int i = 0; i < left.length; i++) {
      final int j = i << 1;
      mArray[j] = left[i];
      mArray[j + 1] = right[i];
    }
  }

  /**
   * Get the underlying value from the array.
   * @param index into the array.
   * @return the value.
   */
  @Override
  public int get(final int index) {
    return mArray[index];
  }

  /**
   * Get the length of the wrapped array.
   * @return the length of the wrapped array.
   */
  @Override
  public int length() {
    return mArray.length;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append(Arrays.toString(mArray));
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mArray != null);
    return true;
  }
}
