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
