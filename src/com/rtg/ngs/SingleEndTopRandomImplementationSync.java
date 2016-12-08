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
 * Synchronized version
 */
public class SingleEndTopRandomImplementationSync extends SingleEndTopRandomImplementation {

  private static final int NUMBER_OF_THREAD_LOCKS = 1 << 16;
  private static final int THREAD_LOCK_MASK = NUMBER_OF_THREAD_LOCKS - 1;

  /* Array of locks for multiple threads */
  private Object[] mThreadLocks;

  /**
   * Constructor
   * @param numReads number of paired reads
   */
  public SingleEndTopRandomImplementationSync(int numReads) {
    super(numReads);
    initThreadLocks();
  }
  SingleEndTopRandomImplementationSync(int numReads, long seed) {
    super(numReads, seed);
    initThreadLocks();
  }

  private void initThreadLocks() {
    mThreadLocks = new Object[NUMBER_OF_THREAD_LOCKS];
    for (int i = 0; i < mThreadLocks.length; ++i) {
      mThreadLocks[i] = new Object();
    }
  }

  @Override
  public void update(int readId, int templateId, int templateStart, boolean reverse, int alignScore, int n) {
    synchronized (mThreadLocks[syncId(readId)]) {
      super.update(readId, templateId, templateStart, reverse, alignScore, n);
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
    return "SingleEndTopRandomImplementationSync";
  }


}
