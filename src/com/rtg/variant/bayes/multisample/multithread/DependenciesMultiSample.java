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
    ORDERING.setLink(JobType.DANGLING, 0, JobType.BED); // Needed to get complexities out, since COMPLEX doesn't return them as a result
    ORDERING.setLink(JobType.COMPLEX, 0, JobType.BED);  // But needs to be after complex calling has put status into the complexities.

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
