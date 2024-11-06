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

package com.rtg.variant.bayes.snp;

import java.util.Arrays;

import com.rtg.variant.bayes.Description;

/**
 * Provides an array-backed set of hypothesis names.
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
    for (int i = 1; i < names.length; ++i) {
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

  @Override
  public String toString() {
    return Arrays.toString(mNames);
  }
}
