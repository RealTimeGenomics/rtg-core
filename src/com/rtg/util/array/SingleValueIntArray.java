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

/**
 * Implementation of integer array that only contains the same value in every position
 */
public class SingleValueIntArray implements ImmutableIntArray {

  private final int mVal;
  private final int mLength;

  /**
   * @param val value to fill array with
   * @param length length of array
   */
  public SingleValueIntArray(int val, int length) {
    mVal = val;
    mLength = length;
  }

  @Override
  public int get(int index) {
    return mVal;
  }

  @Override
  public int length() {
    return mLength;
  }

}
