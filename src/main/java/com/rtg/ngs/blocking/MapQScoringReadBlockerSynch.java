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
 * Synchronized version of MapQScoringReadBlocker.  This is slightly
 * slower than MapQScoringReadBlocker, so should not be used in
 * single-threaded situations.
 *
 */
public class MapQScoringReadBlockerSynch extends MapQScoringReadBlocker {

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
  public MapQScoringReadBlockerSynch(final int count, final int threshold, final String title) {
    super(count, threshold, title);
    mThreadLocks = new Object[NUMBER_OF_THREAD_LOCKS];
    for (int i = 0; i < mThreadLocks.length; ++i) {
      mThreadLocks[i] = new Object();
    }
  }

  /**
   * Creates a counter for <code>count</code> records blocking at <code>
   * threshold</code>.
   *
   * @param count number of reads
   * @param threshold blocking threshold in range 1 to 255.
   */
  public MapQScoringReadBlockerSynch(final int count, final int threshold) {
    this(count, threshold, "multithreaded blocked pairings");
  }

  /**
   * This is not synchronized, so must only be called by one
   * thread, after all other threads have finished using this object.
   */
  @Override
  public void close() {
    super.close();
  }
  @Override
  public int increment(final int r, final int score) {
    // sync the rest of the function
    synchronized (mThreadLocks[r & THREAD_LOCK_MASK]) {
      return super.increment(r, score);
    }
  }

  @Override
  public boolean isBlocked1(final int r, final int score) {
    // sync the rest of the function
    synchronized (mThreadLocks[r & THREAD_LOCK_MASK]) {
      return super.isBlocked1(r, score);
    }
  }

  @Override
  public boolean isBlocked2(final int r, final int score) {
    // sync the rest of the function
    synchronized (mThreadLocks[r & THREAD_LOCK_MASK]) {
      return super.isBlocked2(r, score);
    }
  }

  @Override
  public int getTerminationScore(int r) {
    synchronized (mThreadLocks[r & THREAD_LOCK_MASK]) {
      if (getCount1(r) > 1) {
        return getScore1(r);
      }
    }
    return getScore2(r);
  }

  @Override
  public String toString() {
    return "ScoringReadBlockerSynch";
  }
}
