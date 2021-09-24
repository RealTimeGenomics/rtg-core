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
 * Maps from a {@link JobId} to a {@link Job}. This will vary with each application.
 * @param <J> the type of the job identifiers.
 */
public interface JobFactory<J extends JobId<J>> {

  /**
   * Create the job associated with the identifier.
   * @param id the unique identifier for the job.
   * @param arguments the arguments needed for execution of the job.
   * @return the job.
   */
  Job<J> job(J id, Result[] arguments);
}
