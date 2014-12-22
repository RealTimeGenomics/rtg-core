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

import java.util.HashMap;

/**
 * Utility class for reporting progress from multiple threads using the trailing
 * edge of the reports from individual threads.
 *
 */
public class ParallelProgress {

  private static class MyInteger {
    private int mValue = 0;

    MyInteger(final int v) {
      setValue(v);
    }

    int getValue() {
      return mValue;
    }

    void setValue(int value) {
      mValue = value;
    }

  }

  private final String mName;
  private final HashMap<Thread, MyInteger> mIndividualThreadProgress = new HashMap<>();
  private int mLastReportedProgress = -1;

  /**
   * Create progress reporter.
   * @param name name of progress
   */
  public ParallelProgress(final String name) {
    mName = name;
    Diagnostic.progress("Starting: " + mName);
  }

  /**
   * Update progress safe for multithread.
   * @param value position in current thread
   */
  public synchronized void updateProgress(final int value) {
    final Thread t = Thread.currentThread();
    final MyInteger v = mIndividualThreadProgress.get(t);
    if (v == null) {
      mIndividualThreadProgress.put(t, new MyInteger(value));
    } else {
      v.setValue(value);
    }
    // Only bother to check for global update if this value exceeds last report
    if (value > mLastReportedProgress) {
      int min = Integer.MAX_VALUE;
      for (final MyInteger m : mIndividualThreadProgress.values()) {
        if (m.getValue() < min) {
          min = m.getValue();
        }
      }
      if (min > mLastReportedProgress) {
        mLastReportedProgress = min;
        Diagnostic.progress("Processed " + min + "% of " + mName);
      }
    }
  }

  /**
   * Finish reporting progress.
   */
  public void close() {
    Diagnostic.progress("Finished: " + mName);
  }
}
