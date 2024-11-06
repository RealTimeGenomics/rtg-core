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
