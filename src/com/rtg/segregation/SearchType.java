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
 * What type of link in a chain of blocks.
 */
enum SearchType {
  /** There is a crossover between the previous block and this block. */
  XO("BX", "X"),
  /** This block starts a new sequence disconnected from the previous block. */
  New("BN", "N"),
  /** Treat this block as an error. */
  Error("BE", null),
  /** This block is consistent with the previous block. */
  OK("BL", "N");

  private final String mCode;

  private final String mBedCode;

  /**
   * @param code two letter string used when writing blocks.
   * @param bedCode one letter code used in regions file.
   */
  private SearchType(String code, String bedCode) {
    mCode = code;
    mBedCode = bedCode;
  }

  public String code() {
    return mCode;
  }

  public String bedCode() {
    return mBedCode;
  }
}
