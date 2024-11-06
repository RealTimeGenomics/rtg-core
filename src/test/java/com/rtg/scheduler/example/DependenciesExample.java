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

import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Set;

import com.rtg.scheduler.Dependencies;
import com.rtg.scheduler.LookAhead;
import com.rtg.scheduler.enumtime.EnumTimeId;
import com.rtg.scheduler.enumtime.LocalOrdering;
import com.rtg.util.integrity.IntegralAbstract;

import org.junit.Assert;

/**
 * Example which includes the sort of patterns expected in the multi-sample
 * variant caller.
 */
public class DependenciesExample extends IntegralAbstract implements Dependencies<ExampleJobId> {

  /**
   * Used to compute ordering.
   */
  public static final LocalOrdering<JobType> ORDERING = new LocalOrdering<>(JobType.values());
  static {

    //The order of these calls is important as it determines the order of the arguments to the respective jobs.
    ORDERING.setLink(JobType.INCR, 1, JobType.DANGLING);
    ORDERING.setLink(JobType.INCR, 0, JobType.DANGLING);

    ORDERING.setLink(JobType.BED, 1, JobType.BED);
    ORDERING.setLink(JobType.DANGLING, 0, JobType.BED);

    ORDERING.setLink(JobType.DANGLING, 0, JobType.COMPLEX);

    ORDERING.setLink(JobType.DANGLING, 0, JobType.MERGE);
    ORDERING.setLink(JobType.COMPLEX, 0, JobType.MERGE);

    ORDERING.setLink(JobType.OUT, 1, JobType.OUT);
    ORDERING.setLink(JobType.MERGE, 0, JobType.OUT);

    ORDERING.freeze();
  }

  private final int mNumberChunks;

  private int mNext = 0;

  /**
   * @param numberChunks number of chunks being scheduled.
   */
  public DependenciesExample(final int numberChunks) {
    mNumberChunks = numberChunks;

  }

  @Override
  public int delta() {
    return 1;
  }

  @Override
  public Set<ExampleJobId> from(final ExampleJobId id) {
    final Set<ExampleJobId> res = new LinkedHashSet<>();
    final int chunk = id.time();
    for (final EnumTimeId<JobType> frid : ORDERING.from(id.type())) {
      final int incr = frid.time();
      assert incr <= 0;
      final int fromChunk = chunk + incr;
      res.add(fromChunk < 0 ? null : new ExampleJobId(mNumberChunks, fromChunk, frid.type()));
    }
    return res;
  }

  @Override
  public Collection<ExampleJobId> to(final ExampleJobId id) {
    final Set<ExampleJobId> res = new HashSet<>();
    final int chunk = id.time();
    final Set<EnumTimeId<JobType>> toSet = ORDERING.to(id.type());
    for (final EnumTimeId<JobType> toid : toSet) {
      final int toIncr = toid.time();
      assert toIncr >= 0;
      //System.err.println("toIncr=" + toIncr);
      final int toChunk = chunk + toIncr;
      res.add(toChunk >= mNumberChunks ? null : new ExampleJobId(mNumberChunks, toChunk, toid.type()));
    }
    return res;
  }

  @Override
  public ExampleJobId next(LookAhead lookAhead) {
    final ExampleJobId res;
    if (mNext >= mNumberChunks) {
      res = null;
    } else {
      res = new ExampleJobId(mNumberChunks, mNext, JobType.INCR);
    }
    ++mNext;
    return res;
  }

  @Override
  public boolean integrity() {
    Assert.assertTrue(mNumberChunks > 0);
    Assert.assertTrue(-2 <= mNext && mNext <= mNumberChunks);
    return true;
  }
}
