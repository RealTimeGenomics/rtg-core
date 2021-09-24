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

/**
 * Trivial Executor that does everything sequentially.
 * @param <J> job id
 */
public class ExecutorSequential<J extends JobId<J>> implements Executor<J> {

  private final Scheduler<J> mScheduler;
  /**
   * @param scheduler that ensures no conflicts between jobs.
   */
  public ExecutorSequential(final Scheduler<J> scheduler) {
    assert scheduler.checkEmpty();
    mScheduler = scheduler;
  }

  @Override
  public void run() throws IOException {
    Job<J> current = null;
    Result result = null;
    long time = -1;
    while (true) {
      final Job<J> next = mScheduler.doneNext(current == null ? null : current.id(), result, time);
      if (next == null) {
        break;
      }
      current = next;
      final long t0 = System.nanoTime();
      result = current.runCatch();
      final long t1 = System.nanoTime();
      time = t1 - t0;
    }
    assert mScheduler.checkEmpty();
  }


}
