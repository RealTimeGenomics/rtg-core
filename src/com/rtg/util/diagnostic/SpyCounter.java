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
public class SpyCounter {
  private final String mName;
  private long mCount = 0;

  /**
   * @param name used in reporting results.
   */
  public SpyCounter(String name) {
    mName = name;
    Spy.add(this);
  }

  /**
   * Increment the counter.
   */
  public void increment() {
    mCount++;
  }

  @Override
  public String toString() {
    return mName + " counts " + mCount;
  }
}
