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
 * Synchronize writes to an arbitrary <code>UptoNStore</code>
 */
public class UptoNStoreSync implements UptoNStore {

  private static final int NUMBER_OF_THREAD_LOCKS = 1 << 16;
  private static final int THREAD_LOCK_MASK = NUMBER_OF_THREAD_LOCKS - 1;

  /* Array of locks for multiple threads */
  private final Object[] mThreadLocks;

  private final UptoNStore mEnclosed;

  /**
   * Constructor
   * @param enclosed implementation to enclose
   */
  public UptoNStoreSync(UptoNStore enclosed) {
    mThreadLocks = new Object[NUMBER_OF_THREAD_LOCKS];
    for (int i = 0; i < mThreadLocks.length; i++) {
      mThreadLocks[i] = new Object();
    }
    mEnclosed = enclosed;
  }

  @Override
  public String histogram() {
    return mEnclosed.histogram();
  }

  @Override
  public void setResults(MatchResult results, int encodedReadId) {
    mEnclosed.setResults(results, encodedReadId);
  }

  @Override
  public void process(long templateId, boolean reverse, int encodedReadId, int tStart, int scoreIndel) {
    synchronized (mThreadLocks[syncId(encodedReadId)]) {
      mEnclosed.process(templateId, reverse, encodedReadId, tStart, scoreIndel);
    }
  }

  /**
   * Get the object on which to synchronize.
   * @param encodedReadId read identifier being accessed.
   * @return the index for the read identifier.
   */
  static int syncId(int encodedReadId) {
    return encodedReadId & THREAD_LOCK_MASK;
  }

  @Override
  public String toString() {
    return mEnclosed.toString() + "Sync";
  }


}
