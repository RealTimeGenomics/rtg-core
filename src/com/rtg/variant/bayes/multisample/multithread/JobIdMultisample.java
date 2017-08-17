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

package com.rtg.variant.bayes.multisample.multithread;


import com.rtg.scheduler.JobId;
import com.rtg.scheduler.Result;
import com.rtg.scheduler.enumtime.EnumTimeId;
import com.rtg.util.integrity.Exam;

/**
 */
public class JobIdMultisample extends EnumTimeId<JobType> implements JobId<JobIdMultisample> {

  /**
   * @param numberChunks total number of allowed chunks.
   * @param chunk identifier for the chunk (0 based).
   * @param type of the job.
   */
  public JobIdMultisample(final int numberChunks, final int chunk, final JobType type) {
    super(chunk, type);
    assert chunk >= 0 && chunk <= numberChunks + 1;
  }

  @Override
  public boolean validArguments(Result[] results) {
    return type().validArguments(results);
  }

  @Override
  public boolean validResult(Result result) {
    return type().validResult(result);
  }

  @Override
  public int compareTo(final JobIdMultisample that) {
    final int c0 = this.time() - that.time();
    if (c0 != 0) {
      return c0;
    }
    final boolean before = DependenciesMultiSample.ORDERING.before(this.type(), that.type());
    if (before) {
      return -1;
    }
    final boolean after = DependenciesMultiSample.ORDERING.before(that.type(), this.type());
    if (after) {
      return +1;
    }
    return 0;
  }

  @Override
  public boolean equals(Object arg0) {
    if (!(arg0 instanceof JobIdMultisample)) {
      return false;
    }
    return 0 == compareTo((JobIdMultisample) arg0);
  }

  @Override
  public int hashCode() {
    return super.hashCode();
  }

  @Override
  public boolean integrity() {
    super.integrity();
    Exam.assertTrue(time() >= 0);
    return true;
  }

}
