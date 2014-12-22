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
package com.rtg.reference;

/**
 */
public enum Ploidy {
  /** Is being ignored (because of sex) */
  NONE(0),
  /** Single copy. */
  HAPLOID(1),
  /** Two copies. */
  DIPLOID(2),
  /** Many, non-fixed number of copies. */
  POLYPLOID(-1),
  ;

  private final int mCount;

  private Ploidy(int count) {
    mCount = count;
  }

  /**
   * @return the count of copies, or -1 for indeterminate (i.e. for polyploid).
   */
  public int count() {
    return mCount;
  }
}
