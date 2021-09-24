/*
 * Copyright (c) 2017. Real Time Genomics Limited.
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
 *  Bundles up an R1 / R2 pair
 */
class FastqPair {
  private final FastqSequence mR1;
  private final FastqSequence mR2;

  /**
   * @param r1 the R1 AKA left read
   * @param r2 the R2 AKA right read
   */
  FastqPair(FastqSequence r1, FastqSequence r2) {
    mR1 = r1;
    mR2 = r2;
  }

  /**
   * @return the R1 AKA left read
   */
  public FastqSequence r1() {
    return mR1;
  }

  /**
   * @return the R2 AKA right read
   */
  public FastqSequence r2() {
    return mR2;
  }
}
