/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.ngs;



/**
 * Synchronized version
 */
public class PairedTopRandomImplementationSync extends PairedTopRandomImplementation {

  private static final int NUMBER_OF_THREAD_LOCKS = 1 << 16;
  private static final int THREAD_LOCK_MASK = NUMBER_OF_THREAD_LOCKS - 1;

  /* Array of locks for multiple threads */
  private Object[] mThreadLocks;

  /**
   * Constructor
   * @param numReads number of paired reads
   */
  public PairedTopRandomImplementationSync(int numReads) {
    super(numReads);
    initThreadLocks();
  }
  PairedTopRandomImplementationSync(int numReads, long seed) {
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
  public void update(int readId, int templateId, int templateStart, boolean reverse, int mateTemplateStart, boolean mateReverse, int comboScore, int n) {
    synchronized (mThreadLocks[syncId(readId)]) {
      super.update(readId, templateId, templateStart, reverse, mateTemplateStart, mateReverse, comboScore, n);
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
    return "PairedTopRandomImplementationSync";
  }


}
