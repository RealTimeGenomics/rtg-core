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
package com.rtg.variant.realign;

/**
 * Accrue counts for homopolymer regions
 */
public class HomopolymerRepeats {

  private final int[] mForward;
  private final int[] mReverse;

  /**
   * @param sequence bytes for template using ordinal for values
   */
  public HomopolymerRepeats(ByteArrayAdaptor sequence) {
    this(sequence, 0, sequence.length());
  }

  /**
   * @param sequence bytes (0 = N, ... 4 = T).
   * @param start position on the sequence (0 based).
   * @param end position on sequence (0 based, exclusive).
   */
  public HomopolymerRepeats(ByteArrayAdaptor sequence, final int start, final int end) {
    final int length = end - start;
    mForward = new int[length];
    mReverse = new int[length];
    int startRepeat = 0;
    byte last = -1; //initial value doesn't matter
    int count = 0;
    for (int i = start, j = 0; i < end; i++, j++) {
      if (j != startRepeat) {
        if (last != sequence.get(i)) {
          mForward[j - 1] = count;
          mReverse[startRepeat] = count;
          count = 0;
          startRepeat = j;
        }
        //everything starts at 0 in java so don't need to zero in else
      }
      count++;
      last = sequence.get(i);
    }
    mForward[mForward.length - 1] = count;
    mReverse[startRepeat] = count;
  }

  /**
   * A forward count. A position with a value of zero represents inside a homopolymer region
   * any other number specifies the length of the homopolymer region that just finished (the last nucleotide of
   * the repeat is at <code>i</code>).
   * @param i the position to get the value (0 based).
   * @return the count.
   */
  public int forward(int i) {
    if (i < 0 || i >= mForward.length) {
      return -1;
    }
    return mForward[i];
  }

  /**
   * A reverse count. A position with a value of zero represents inside a homopolymer region
   * any other number specifies the length of the homopolymer region that just finished (the last nucleotide of
   * the repeat is at <code>i</code>).
   * @param i the position to get the value (0 based).
   * @return the count.
   */
  public int reverse(int i) {
    if (i < 0 || i >= mReverse.length) {
      return -1;
    }
    return mReverse[i];
  }
}
