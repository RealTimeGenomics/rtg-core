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
