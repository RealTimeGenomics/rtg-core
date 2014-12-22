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

package com.rtg.variant.bayes.snp;

import com.rtg.variant.bayes.Description;

/**
 * Provides an array-backed set of hypothesis names
 */
public class DescriptionCommon extends Description {

  private final String[] mNames;

  private final int mMinLength;

  private final int mMaxLength;

  /**
   * Make one.
   * @param names for each hypothesis.
   */
  public DescriptionCommon(final String... names) {
    if (names.length == 0) {
      throw new IllegalArgumentException();
    }
    mNames = names;
    int min = names[0].length();
    int max = names[0].length();
    for (int i = 1; i < names.length; i++) {
      final int len = names[i].length();
      if (len < min) {
        min = len;
      }
      if (len > max) {
        max = len;
      }
    }
    mMinLength = min;
    mMaxLength = max;
  }

  @Override
  public String name(int index) {
    return mNames[index];
  }

  @Override
  public int size() {
    return mNames.length;
  }

  @Override
  public int minLength() {
    return mMinLength;
  }

  @Override
  public int maxLength() {
    return mMaxLength;
  }

}
