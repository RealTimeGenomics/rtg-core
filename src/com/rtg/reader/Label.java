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

package com.rtg.reader;

/**
 * Holds sequence names as two parts, an initial short name and a suffix.
 */
public class Label {
  private final String mLabel;
  private final String mSuffix;

  /**
   * Creates a label loading the fields.
   *
   * @param label  the short name for the label
   * @param suffix the suffix for the label
   */
  public Label(String label, String suffix) {
    mLabel = label == null ? "" : label;
    mSuffix = suffix == null ? "" : suffix;
  }

  /**
   * Returns the short name for the label.
   *
   * @return short name
   */
  public String label() {
    return mLabel;
  }

  /**
   * Returns the suffix for the label.
   *
   * @return suffix
   */
  public String suffix() {
    return mSuffix;
  }

  /**
   * Returns the full name for the label.
   *
   * @return full name
   */
  public String fullLabel() {
    return mLabel + mSuffix;
  }

  /**
   * Returns length of label.
   *
   * @return length of label
   */
  public int length() {
    return mLabel.length() + mSuffix.length();
  }

  @Override
  public String toString() {
    return fullLabel();
  }
}
