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

import com.rtg.scheduler.JobId;
import com.rtg.scheduler.Result;
import com.rtg.scheduler.enumtime.EnumTimeId;

import org.junit.Assert;

/**
 */
public class ExampleJobId extends EnumTimeId<JobType>  implements JobId<ExampleJobId> {

  /**
   * @param numberChunks total number of allowed chunks.
   * @param chunk identifier for the chunk (0 based).
   * @param type of the job.
   */
  public ExampleJobId(final int numberChunks, final int chunk, final JobType type) {
    super(chunk, type);
    assert chunk >= 0 && chunk < numberChunks;
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
  public int compareTo(final ExampleJobId that) {
    final int c0 = this.time() - that.time();
    if (c0 != 0) {
      return c0;
    }
    final boolean before = DependenciesExample.ORDERING.before(this.type(), that.type());
    if (before) {
      return -1;
    }
    final boolean after = DependenciesExample.ORDERING.before(that.type(), this.type());
    if (after) {
      return +1;
    }
    return 0;
  }

  @Override
  public boolean equals(Object arg0) {
    return super.equals(arg0);
  }

  @Override
  public int hashCode() {
    return super.hashCode();
  }

  @Override
  public boolean integrity() {
    super.integrity();
    Assert.assertTrue(time() >= 0);
    return true;
  }

}
