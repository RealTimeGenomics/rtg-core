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
package com.rtg.sam;

/**
 * Records consensus statistics
 */
public class PileUp {

  private final int mTemplateLength;
  private long mCount = 0;
  private long mNCount = 0;

  private final int[] mA;
  private final int[] mG;
  private final int[] mT;
  private final int[] mC;

  PileUp(final int templateLength) {
    mTemplateLength = templateLength;
    mA = new int[templateLength];
    mG = new int[templateLength];
    mT = new int[templateLength];
    mC = new int[templateLength];
  }

  void add(char base, int position) {
    switch (Character.toLowerCase(base)) {
    case 'a':
      mA[position]++;
      break;
    case 'g':
      mG[position]++;
      break;
    case 't':
      mT[position]++;
      break;
    case 'c':
      mC[position]++;
      break;
    case 'n':
      ++mNCount;
      break;
    default:
      break;
    }
    ++mCount;
  }

  long consensus() {
    long c = mNCount;
    for (int i = 0; i < mTemplateLength; ++i) {
      c += Math.max(mA[i], Math.max(mG[i], Math.max(mT[i], mC[i])));
    }

    return c;
  }

  long total() {
    return mCount;
  }

  double coverage() {
    return total() / (double) mTemplateLength;
  }
}
