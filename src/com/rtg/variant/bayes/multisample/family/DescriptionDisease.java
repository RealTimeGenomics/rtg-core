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

package com.rtg.variant.bayes.multisample.family;

import com.rtg.variant.bayes.Description;

/**
 */
public class DescriptionDisease extends Description {

  private static final String NONE = "NONE";

  private final Description mDescription;

  /**
   * @param description underlying description to which NONE is prepended.
   */
  DescriptionDisease(Description description) {
    super();
    mDescription = description;
  }

  @Override
  public String name(int index) {
    if (index == 0) {
      return NONE;
    }
    return mDescription.name(index - 1);
  }

  @Override
  public int size() {
    return mDescription.size() + 1;
  }

  @Override
  public int minLength() {
    return mDescription.minLength();
  }

  @Override
  public int maxLength() {
    return mDescription.maxLength();
  }
}
