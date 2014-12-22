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

package com.rtg.util.diagnostic;

/**
 * Collect counts and have them reported when <code>Spy.report()</code> is called.
 */
public class SpyHistogram {
  private final String mName;
  private long mCount = 0;
  private final long[] mHisto;

  /**
   * @param name used in reporting results.
   * @param length number of counters in histogram.
   */
  public SpyHistogram(final String name, final int length) {
    mName = name;
    mHisto = new long[length];
    Spy.add(this);
  }

  /**
   * Increment the histogram counter as determined by index.
   * If beyond the specified length then increment a separate counter.
   * @param index specifies counter to be incremented assumed to be &ge; 0.
   */
  public void increment(final int index) {
    if (index >= mHisto.length) {
      mCount++;
      return;
    }
    mHisto[index]++;
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append(mName).append(" [").append(mHisto.length).append("] ");
    for (long aMHisto : mHisto) {
      sb.append(aMHisto).append(" ");
    }
    sb.append("...").append(mCount);
    return sb.toString();
  }
}
