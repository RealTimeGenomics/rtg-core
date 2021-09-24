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
package com.rtg.ngs.blocking;

/**
 * Synchronized version of <code>ReadBlocker</code>
 */
public class ReadBlockerSync extends ReadBlocker {
  private static final int NUMBER_OF_THREAD_LOCKS = 1 << 16;
  private static final int THREAD_LOCK_MASK = NUMBER_OF_THREAD_LOCKS - 1;

  /* Array of locks for multiple threads */
  private final Object[] mThreadLocks;

  /**
   * Creates a counter for <code>count</code> records blocking at <code>
   * threshold</code>.
   *
   * @param count number of reads
   * @param threshold blocking threshold in range 1 to 255
   * @param title a title to use during logging
   */
  public ReadBlockerSync(final long count, final int threshold, final String title) {
    super(count, threshold, title);
    mThreadLocks = new Object[NUMBER_OF_THREAD_LOCKS];
    for (int i = 0; i < mThreadLocks.length; ++i) {
      mThreadLocks[i] = new Object();
    }
  }
  @Override
  public void increment(final int r) {
    // sync the rest of the function
    synchronized (mThreadLocks[r & THREAD_LOCK_MASK]) {
      if (!super.isBlocked(r)) {
        super.increment(r);
      }
    }
  }

  @Override
  public boolean isBlocked(final int r) {
    // sync the rest of the function
    synchronized (mThreadLocks[r & THREAD_LOCK_MASK]) {
      return super.isBlocked(r);
    }
  }

}
