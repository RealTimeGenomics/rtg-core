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
import com.rtg.scheduler.JobFactory;
import com.rtg.scheduler.Result;
import com.rtg.util.integrity.IntegralAbstract;

import junit.framework.Assert;

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
