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

package com.rtg.reader;

/**
 * Clips sequences based on average quality in a sliding window.
 *
 */
public class GilletteReadTrimmer implements ReadTrimmer {

  private final int mWindowSize;
  private final int mQualityThreshold;

  /**
   * Construct a Gillette read trimmer
   * @param windowSize size of the window to look over
   * @param qualityThreshold the threshold the window must average higher than
   */
  public GilletteReadTrimmer(int windowSize, int qualityThreshold) {
    mWindowSize = windowSize;
    mQualityThreshold = qualityThreshold;
  }

  @Override
  public int getTrimPosition(byte[] qualities, int length) {
    final int[] quals = new int[mWindowSize];
    int cutoffIndex = length;
    double sum = 0.0;
    for (int i = 0; i < cutoffIndex; i++) {
      if (i >= quals.length) {
        if (sum / quals.length < mQualityThreshold) {
          cutoffIndex = i;
        }
      }
      final int i2 = i % quals.length;
      sum -= quals[i2];
      quals[i2] = (int) qualities[i];
      sum += quals[i2];
    }
    if (cutoffIndex > mWindowSize) {
      return cutoffIndex;
    }
    return 0;
  }

}
