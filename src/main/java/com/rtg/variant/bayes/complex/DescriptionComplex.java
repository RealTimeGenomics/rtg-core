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

package com.rtg.variant.bayes.complex;

import java.util.List;

import com.rtg.variant.bayes.Description;
import com.rtg.variant.match.Match;

/**
 * Class containing static description for all the hypothesis
 */
public class DescriptionComplex extends Description {

  private final List<? extends Match> mNames;

  private final int mMinLength;

  private final int mMaxLength;

  /**
   * Make one.
   * @param matches corresponding to each hypothesis.
   */
  public DescriptionComplex(List<? extends Match> matches) {
    mNames = matches;
    int min = Integer.MAX_VALUE;
    int max = 0;
    for (final Match match : matches) {
      final int length = match.length();
      min = Math.min(min, length);
      max = Math.max(max, length);
    }
    mMinLength = min;
    mMaxLength = max;
  }

  /**
   * Return match at given index
   * @param index index of match
   * @return a match
   */
  public Match match(int index) {
    return mNames.get(index);
  }

  @Override
  public String name(int index) {
    return mNames.get(index).readString();
  }

  @Override
  public int size() {
    return mNames.size();
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
    return mNames.toString();
  }
}
