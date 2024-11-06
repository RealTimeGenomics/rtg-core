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
import com.rtg.scheduler.JobFactory;
import com.rtg.scheduler.Result;
import com.rtg.util.integrity.IntegralAbstract;

import org.junit.Assert;

/**
 * Example which includes the sort of patterns expected in the multi-sample
 * variant caller.
 */
public class JobFactoryExample extends IntegralAbstract implements JobFactory<ExampleJobId> {

  /** For testing result in bed output. */
  private final StringBuilder mBed = new StringBuilder();

  /** For testing result in out output. */
  private final StringBuilder mOut = new StringBuilder();

  private final int mNumberChunks;

  /**
   * @param numberChunks number chunks.
   */
  public JobFactoryExample(final int numberChunks) {
    super();
    mNumberChunks = numberChunks;
  }

  /**
   * @return For testing result in bed output.
   */
  public StringBuilder bed() {
    return mBed;
  }

  /**
   * @return For testing result in out output.
   */
  public StringBuilder out() {
    return mOut;
  }

  @Override
  public Job<ExampleJobId> job(final ExampleJobId id, final Result[] arguments) {
    final StringBuilder res;
    if (id.time() == mNumberChunks - 1) {
      if (id.type() == JobType.BED) {
        res = mBed;
      } else if (id.type() == JobType.OUT) {
        res = mOut;
      } else {
        res = null;
      }
    } else {
      res = null;
    }
    return new MockJob(id, arguments, res);
  }

  @Override
  public boolean integrity() {
    Assert.assertTrue(mNumberChunks > 0);
    return true;
  }
}
