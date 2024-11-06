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

package com.rtg.variant.bayes.multisample.multithread;

import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import com.rtg.scheduler.Dependencies;
import com.rtg.scheduler.LookAhead;
import com.rtg.scheduler.enumtime.EnumTimeId;
import com.rtg.scheduler.enumtime.LocalOrdering;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Example which includes the sort of patterns expected in the multi-sample
 * variant caller.
 */
public class DependenciesMultiSample extends IntegralAbstract implements Dependencies<JobIdMultisample> {

  /**
   * Used to compute ordering.
   */
  public static final LocalOrdering<JobType> ORDERING = new LocalOrdering<>(JobType.values());
  static {

    // This helps to visualize the dependencies. But also, see the special cases implemented in from() and to().

    // dot -Tpdf -o jobdeps.pdf <(grep set[L]ink src/com/rtg/variant/bayes/multisample/multithread/DependenciesMultiSample.java | grep -v "^ \+//" | sed -e "s|.*(\(.*\)).*|\1|g" -e "s|JobType.||g" | tr -d ',' | awk 'BEGIN{print "digraph { \nlabelloc=\"t\"\nlabel=\"Scheduler Job Dependencies\""}{print "    "$1"_0 -> "$3"_"$2" [label=\""$2"\"]"; nodes[$1"_0"]=1; nodes[$3"_"$2]=1; if ($2=="0") print "    "$1"_1 -> "$3"_1 [label=\"0\"]"} END{ for (node in nodes) print "    "node" [label=\""node"\"]"; print "}"}')


    //The order of these calls is important as it determines the order of the arguments to the respective jobs.
    ORDERING.setLink(JobType.INCR, 1, JobType.DANGLING);
    ORDERING.setLink(JobType.INCR, 0, JobType.DANGLING);
    ORDERING.setLink(JobType.DANGLING, 1, JobType.DANGLING);  // Dangling adjustment propagation through consecutive chunks

    ORDERING.setLink(JobType.BED, 1, JobType.BED);      // Ensure output ordering
    ORDERING.setLink(JobType.COMPLEX, 0, JobType.BED);

    ORDERING.setLink(JobType.DANGLING, 0, JobType.COMPLEX); // DANGLING passes the region from the previous timestep to complex

    ORDERING.setLink(JobType.INCR, 1, JobType.FLUSH); // Pass through the complexities so we know what region to flush.
    ORDERING.setLink(JobType.COMPLEX, 0, JobType.FLUSH); // Ensure DANGLING has adjusted region and COMPLEX has finished calling

    ORDERING.setLink(JobType.INCR, 1, JobType.FILTER);   // Only needed to propagate max read length for equivalent filtering limiting
    ORDERING.setLink(JobType.COMPLEX, 0, JobType.FILTER);

    ORDERING.setLink(JobType.OUT, 1, JobType.OUT); // Ensure output ordering and needed for equivalent filtering across boundaries
    ORDERING.setLink(JobType.FILTER, 0, JobType.OUT);
    ORDERING.setLink(JobType.FLUSH, 0, JobType.OUT);

    ORDERING.freeze();
  }

  private final int mNumberChunks;

  private int mNext = 0;

  /**
   * @param numberChunks number chunks.
   */
  public DependenciesMultiSample(final int numberChunks) {
    mNumberChunks = numberChunks;
  }


  @Override
  public int delta() {
    //Maximum dependency depth for multi-sample variant caller is 1 chunk away
    return 1;
  }

  @Override
  public Collection<JobIdMultisample> from(final JobIdMultisample id) {
    final List<JobIdMultisample> res = new LinkedList<>();
    final int chunk = id.time(); //absolute time
    final JobType type = id.type();
    final Set<EnumTimeId<JobType>> fromId = ORDERING.from(id.type());
    //System.err.println(fromId);
    for (final EnumTimeId<JobType> frid : fromId) {
      final int incr = frid.time(); //relative time
      assert incr <= 0;
      final int fromChunk = chunk + incr;
      final JobIdMultisample from;
      final JobType fridType = frid.type();
      if (fromChunk < 0 || ((type == JobType.DANGLING || type == JobType.FLUSH) && chunk == mNumberChunks && fridType == JobType.INCR && incr == 0)) {
        from = null;
      } else {
        final boolean none;
        if (chunk == mNumberChunks + 1) {
          switch (type) {
            case FILTER:
              none = false;
              break;
            case OUT:
              none = fridType == JobType.FILTER || fridType == JobType.FLUSH;
              break;
            default:
              throw new RuntimeException(id.toString());
          }
        } else {
          none = false;
        }
        from = none ? null : new JobIdMultisample(mNumberChunks, fromChunk, fridType);
      }
      res.add(from);
    }
    //System.err.println("from id:" + id + " <- " + res);
    return res;
  }

  @Override
  public Collection<JobIdMultisample> to(final JobIdMultisample id) {
    final Set<JobIdMultisample> res = new HashSet<>();
    final int chunk = id.time();
    final Set<EnumTimeId<JobType>> toSet = ORDERING.to(id.type());
    for (final EnumTimeId<JobType> toid : toSet) {
      final int toIncr = toid.time();
      assert toIncr >= 0;
      final int toChunk = chunk + toIncr;
      final JobType toType = toid.type();
      final boolean inRange;
      if (toChunk <= mNumberChunks) {
        inRange = true;
      } else {
        inRange = toChunk == mNumberChunks + 1 && (toType == JobType.OUT);
      }
      if (inRange) {
        res.add(new JobIdMultisample(mNumberChunks, toChunk, toType));
      }
    }
    return res;
  }

  @Override
  public JobIdMultisample next(LookAhead lookAhead) {
    final JobIdMultisample res;
    if (mNext >= mNumberChunks || !lookAhead.ok(mNext, 0)) {
      res = null;
    } else {
      res = new JobIdMultisample(mNumberChunks, mNext, JobType.INCR);
      ++mNext;
    }
    return res;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mNumberChunks > 0);
    Exam.assertTrue(0 <= mNext && mNext <= mNumberChunks);
    return true;
  }

}
