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
import java.util.concurrent.atomic.AtomicInteger;

import com.rtg.util.IORunnable;
import com.rtg.util.ProgramState;
import com.rtg.util.SimpleThreadPool;

/**
 * Trivial Executor that does everything sequentially.
 * @param <J> job id
 */
public class ExecutorThreaded<J extends JobId<J>> implements Executor<J> {
  private final SimpleThreadPool mPool;
  private final AtomicInteger mThreadsRunning = new AtomicInteger();
  private final Scheduler<J> mScheduler;
  private final int mN;

  /**
   * @param scheduler that ensures no conflicts between jobs.
   * @param n number of pseudo-parallel threads to run.
   */
  public ExecutorThreaded(final Scheduler<J> scheduler, final int n) {
    mPool = new SimpleThreadPool(n, "ExecutorThreaded", false);
    mScheduler = scheduler;
    assert mScheduler.checkEmpty();
    mN = n;
    mThreadsRunning.set(0);
  }

  /**
   * @throws IOException whenever.
   */
  @Override
  public void run() throws IOException {
    //start one thread with null job
    mPool.execute(new ExecutorThread(null));
    //System.err.println("Thread pool terminate");
    mPool.terminate();
    assert mThreadsRunning.get() == 0;
    assert mScheduler.checkEmpty();
  }

  class ExecutorThread implements IORunnable {
    private final Job<J> mStartJob; //can be null

    /**
     * @param startJob initial job to start running.
     */
    ExecutorThread(Job<J> startJob) {
      mStartJob = startJob;
    }

    @Override
    public void run() throws IOException {
      //System.err.println("ExecutorThreaded run()");
      mThreadsRunning.incrementAndGet();
      Job<J> current = mStartJob;
      long nanoTime = -1;
      Result result = null;
      if (current != null) {
        final long t0 = System.nanoTime();
        result = current.runCatch();
        final long t1 = System.nanoTime();
        nanoTime = t1 - t0;
      }
      while (true) {
        //before we start a new job check if there are exceptions from other threads
        ProgramState.checkAbort();
        final Job<J> next = mScheduler.doneNext(current == null ? null : current.id(), result, nanoTime);
        if (next == null) {
          final int running = mThreadsRunning.decrementAndGet();
          assert running >= 0;
          break;
        }
        if (mThreadsRunning.get() < mN) {
          final Job<J> alternate = mScheduler.doneNext(null, null, -1);
          if (alternate != null) {
            mPool.execute(new ExecutorThread(alternate));
          }
        }
        current = next;
        final long t0 = System.nanoTime();
        result = current.runCatch();
        final long t1 = System.nanoTime();
        nanoTime = t1 - t0;
      }
    }
  }
}
