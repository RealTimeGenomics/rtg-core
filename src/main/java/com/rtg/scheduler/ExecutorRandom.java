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

package com.rtg.scheduler;

import java.io.IOException;
import java.util.Random;

/**
 * Trivial Executor that does everything sequentially.
 * @param <J> job id
 */
public class ExecutorRandom<J extends JobId<J>> implements Executor<J> {
  private final Scheduler<J> mScheduler;
  private final int mN;
  private final Random mRandom;
  private final Job<J>[] mRunning;
  private final Result[] mResults;
  private final long[] mTimes;
  private int mNextThread;

  /**
   * @param scheduler that ensures no conflicts between jobs.
   * @param n number of pseudo-parallel threads to run.
   * @param seed for random number generator.
   */
  public ExecutorRandom(final Scheduler<J> scheduler, final int n, final Long seed) {
    assert scheduler.checkEmpty();
    mScheduler = scheduler;
    mN = n;
    mRandom = seed == null ? new Random() : new Random(seed);
    @SuppressWarnings("unchecked")
    final Job<J>[] jobs = (Job<J>[]) new Job<?>[n];
    mRunning = jobs;
    mResults = new Result[n];
    mTimes = new long[n];
    mNextThread = 0;
  }

  @Override
  public void run() throws IOException {
    //true iff the id and results at mNextThread are to be sent
    boolean send = false;
    while (true) {
      final Job<J> next = !send ? mScheduler.doneNext(null, null, -1) : mScheduler.doneNext(mRunning[mNextThread].id(), mResults[mNextThread], mTimes[mNextThread]);
      if (next == null && mNextThread == 0) {
        break;
      }
      if (next == null) {
        send = true;
        choose();
      } else {
        final long t0 = System.nanoTime();
        final Result res = next.runCatch();
        final long t1 = System.nanoTime();
        mTimes[mNextThread] = t1 - t0;
        mRunning[mNextThread] = next;
        mResults[mNextThread] = res;
        ++mNextThread;
        if (mNextThread == mN) {
          send = true;
          choose();
        } else {
          send = false;
        }
      }
    }
    assert mScheduler.checkEmpty();
  }

  /**
   * Choose one item at random to be sent and swap it with the last
   * item at <code>mNextThread</code> which is decremented by 1.
   */
  private void choose() {
    assert mNextThread >= 1;
    final int choose = mRandom.nextInt(mNextThread);
    --mNextThread;
    final Job<J> tempJob = mRunning[mNextThread];
    final Result tempRes = mResults[mNextThread];
    mRunning[mNextThread] = mRunning[choose];
    mResults[mNextThread] = mResults[choose];
    mRunning[choose] = tempJob;
    mResults[choose] = tempRes;
  }
}
