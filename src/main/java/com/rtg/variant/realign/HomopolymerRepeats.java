/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    for (int i = start, j = 0; i < end; ++i, ++j) {
      if (j != startRepeat) {
        if (last != sequence.get(i)) {
          mForward[j - 1] = count;
          mReverse[startRepeat] = count;
          count = 0;
          startRepeat = j;
        }
        //everything starts at 0 in java so don't need to zero in else
      }
      ++count;
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
