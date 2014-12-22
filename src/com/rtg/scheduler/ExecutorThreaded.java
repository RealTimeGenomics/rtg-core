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
    public ExecutorThread(Job<J> startJob) {
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
