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
 * Accumulate timing statistics for scheduled jobs.
 * @param <J> the type of the job identifiers.
 */
public interface JobStatistics<J extends JobId<? extends JobId<?>>> {

  /**
   * Record one execution of id.
   * @param id identifier of the job.
   * @param nanoTime time in nanoseconds it took to run.
   */
  void increment(J id, long nanoTime);

}
