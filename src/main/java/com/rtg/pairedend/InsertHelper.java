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
