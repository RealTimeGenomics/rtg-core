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


/**
 * Takes requests for <code>Job</code>s that can be executed in parallel.
 * @param <J> the type of the job identifiers.
 */
public interface Scheduler<J extends JobId<? extends JobId<?>>> {

  /**
   * Called from an executor to report the completion of a job and to request
   * another job to run. May be called with id null, which indicates a request
   * for a new additional job. This will happens number of times at the start and possibly
   * later. The return result may be null which indicates that no jobs are currently
   * available for execution. If the executor is not currently running any jobs
   * and a null is returned then the executor should terminate.
   * @param id the unique identifier of the job that has just completed (or if null a request for additional work).
   * @param result the result of executing the job uniquely identified by <code>id</code>.
   * @param nanoTime number of nanoseconds taken to run the job.
   * @return a new job to be executed (or null if none are available).
   */
  Job<J> doneNext(J id, Result result, long nanoTime);

  /**
   * Check that everything is empty (used in assertions at start and end).
   * @return true so can be called from an assertion.
   */
  boolean checkEmpty();

  /**
   * Get the look ahead from the scheduler.
   * @return the look ahead object used by the scheduler.
   */
  LookAhead lookAhead();
}
