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

package com.rtg.segregation;

/**
 * Helper for printing set style output.
 */
public class InnerAppend {

  /**
   * @return suitable for outputting a set.
   */
  public static InnerAppend innerSet() {
    return new InnerAppend("{", ", ", "}");
  }

  private final String mPrefix;
  private final String mSeparator;
  private final String mSuffix;

  private boolean mFirstDone = false;

  /**
   * @param prefix done once at start.
   * @param separator before each appended item except the first.
   * @param suffix done once at end.
   */
  public InnerAppend(String prefix, String separator, String suffix) {
    mPrefix = prefix;
    mSeparator = separator;
    mSuffix = suffix;
  }

  /**
   * @return the prefix.
   */
  public String start() {
    assert !mFirstDone;
    return mPrefix;
  }

  /**
   * @param str to be returned with any necessary prepended constant.
   * @return str with prepended constant as necessary.
   */
  public String inner(final String str) {
    if (mFirstDone) {
      return mSeparator + str;
    }
    mFirstDone = true;
    return str;
  }

  /**
   * @return the suffix.
   */
  public String end() {
    return mSuffix;
  }
}
