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
