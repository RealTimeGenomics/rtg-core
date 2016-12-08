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
package com.rtg.ngs;


/**
 * Keeps track of information for determination of unmapped etymology
 */
public class ReadStatusTrackerSync extends ReadStatusTracker {

  //Sync locks for status array in superclass object
  private static final int NUMBER_OF_THREAD_LOCKS = 1 << 16;
  private static final int THREAD_LOCK_MASK = NUMBER_OF_THREAD_LOCKS - 1;

  /* Array of locks for multiple threads */
  private final Object[] mThreadLocks;

  /**
   * @param numReads the number of reads
   * @param stats the statistics tracking object
   */
  public ReadStatusTrackerSync(int numReads, MapStatistics stats) {
    super(numReads, stats);
    mThreadLocks = new Object[NUMBER_OF_THREAD_LOCKS];
    for (int i = 0; i < mThreadLocks.length; ++i) {
      mThreadLocks[i] = new Object();
    }
  }

  /**
   * Add status to read, used when syncing needed
   * @param readId read id
   * @param attr attribute to set
   */
  @Override
  public void addStatus(int readId, int attr) {
    if (isSet(mReadIdStatus[readId], attr)) {
      return;
    }
    synchronized (mThreadLocks[readId & THREAD_LOCK_MASK]) {
      mReadIdStatus[readId] |= attr;
    }
  }
}
