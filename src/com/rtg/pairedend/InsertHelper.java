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
package com.rtg.pairedend;

/**
 * This class calculates the insert sizes for given pair.
 */
public final class InsertHelper {

  private InsertHelper() { }

  /**
   * Calculate the fragment length for a given pair.
   * @param templateStart1 template start for first of the pair
   * @param readLen1 read length for first
   * @param templateStart2 template start for second of the pair
   * @param readLen2 read length for second
   * @return fragment length
   */
  public static int calculateFragmentLength(int templateStart1, int readLen1, int templateStart2, int readLen2) {
    if (templateStart1 < templateStart2) {
      return calculateCorrectFragmentLength(templateStart1, readLen1, templateStart2, readLen2);
    } else {
      return calculateCorrectFragmentLength(templateStart2, readLen2, templateStart1, readLen1);
    }
  }

  /**
   * Calculate the fragment length for a given pair.
   * @param first1 true iff the read referred to in <code>templateStart1</code> and <code>readLen1</code> is first of pair.
   * @param templateStart1 template start for first of the pair
   * @param readLen1 read length for first
   * @param templateStart2 template start for second of the pair
   * @param readLen2 read length for second
   * @return fragment length
   */
  public static int tlen(boolean first1, int templateStart1, int readLen1, int templateStart2, int readLen2) {
    final boolean first = templateStart1 < templateStart2 || (templateStart1 == templateStart2 && first1);
    if (first) {
      return calculateCorrectFragmentLength(templateStart1, readLen1, templateStart2, readLen2);
    } else {
      return -calculateCorrectFragmentLength(templateStart2, readLen2, templateStart1, readLen1);
    }
  }

  private static int calculateCorrectFragmentLength(int leftStart, int leftLength, int rightStart, int rightLength) {
    //This takes into account a situation where the leftmost alignment fully envelops the rightmost alignment
    assert leftStart <= rightStart && leftLength >= 0 && rightLength >= 0;
    return Math.max(leftLength, (rightStart + rightLength) - leftStart);
  }
}
