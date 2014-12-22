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

import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Set;

import com.rtg.scheduler.Dependencies;
import com.rtg.scheduler.LookAhead;
import com.rtg.scheduler.enumtime.EnumTimeId;
import com.rtg.scheduler.enumtime.LocalOrdering;
import com.rtg.util.integrity.IntegralAbstract;

import junit.framework.Assert;

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
    mNext++;
    return res;
  }

  @Override
  public boolean integrity() {
    Assert.assertTrue(mNumberChunks > 0);
    Assert.assertTrue(-2 <= mNext && mNext <= mNumberChunks);
    return true;
  }
}
