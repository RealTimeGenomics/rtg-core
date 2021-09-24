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

import java.util.Arrays;

import com.rtg.scheduler.JobStatistics;
import com.rtg.util.diagnostic.SpyTimer;

/**
 */
public class MultisampleStatistics implements JobStatistics<JobIdMultisample> {

  private final SpyTimer[] mSpies = new SpyTimer[JobType.values().length];
  {
    for (final JobType type : JobType.values()) {
      mSpies[type.ordinal()] = new SpyTimer(type.toString());
    }
  }

  @Override
  public void increment(JobIdMultisample id, long nanoTime) {
    mSpies[id.type().ordinal()].increment(nanoTime);
  }

  @Override
  public String toString() {
    return Arrays.toString(mSpies);
  }

}
