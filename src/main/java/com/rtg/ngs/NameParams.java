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
package com.rtg.ngs;

import com.reeltwo.jumble.annotations.TestClass;

/**
 * Small parameter class for keeping track of need to load names from sequence files
 */
@TestClass("com.rtg.ngs.MapFCliTest")
public class NameParams {
  private final boolean mIncludeNames;
  private final boolean mIncludeFullNames;

  /**
   * @param includeNames true to load names from disk
   * @param includeFullNames true to load full names from disk
   */
  public NameParams(boolean includeNames, boolean includeFullNames) {
    mIncludeNames = includeNames;
    mIncludeFullNames = includeFullNames;
  }

  /**
   * @return true to load names from disk
   */
  public boolean includeNames() {
    return mIncludeNames;
  }

  /**
   * @return true to load full names from disk
   */
  public boolean includeFullNames() {
    return mIncludeFullNames;
  }

}
