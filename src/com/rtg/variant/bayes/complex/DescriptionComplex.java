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
}
