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

package com.rtg.scheduler.example;

import com.rtg.scheduler.Job;
import com.rtg.scheduler.Result;
import com.rtg.util.integrity.Exam;

class MockJob extends Job<ExampleJobId> {

  private final Result[] mArguments;

  private final StringBuilder mResult;

  /**
   * @param jobId the inique identifier for this job.
   * @param arguments the arguments supplied by other jobs.
   * @param result put the result of the job in here as a side effect (demonstrates external communication with the world outside the jobs).
   */
  MockJob(final ExampleJobId jobId, final Result[] arguments, final StringBuilder result) {
    super(jobId);
    mArguments = arguments;
    mResult = result;
  }

  @Override
  public Result run() {
    final StringBuilder sb = new StringBuilder();
    sb.append(id().toString()).append("(");
    for (int i = 0; i < mArguments.length; ++i) {
      if (i > 0) {
        sb.append(", ");
      }
      sb.append(mArguments[i]);
    }
    sb.append(")");
    final String str = sb.toString();
    if (mResult != null) {
      mResult.append(str);
      return null;
    }
    return new MockResult(str);
  }

  @Override
  public boolean integrity() {
    super.integrity();
    Exam.assertNotNull(mArguments);
    return true;
  }

  @Override
  public String toString() {
    return run().toString();
  }

}
