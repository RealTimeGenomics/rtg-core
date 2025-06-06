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

package com.rtg.scheduler;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import com.rtg.util.MathUtils;
import com.rtg.util.Utils;
import com.rtg.util.array.ArrayUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.util.integrity.Integrity;
import com.rtg.variant.bayes.multisample.multithread.JobIdMultisample;
import com.rtg.variant.bayes.multisample.multithread.JobType;

/**
 * A scheduler for jobs that includes a facility to trace the execution.
 * It is synchronized and so will provide a single point bottleneck for
 * any executor using it. This will be less efficient as the number of threads increases
 * or the size of the jobs decreases. (An implementation that avoids this problem appears to be much
 * trickier to implement).
 * @param <J> the type of the job identifiers.
 */
public class SchedulerSynchronized<J extends JobId<J>> implements Scheduler<J>, Integrity {

  private static class CountResult extends IntegralAbstract {
    private int mRefCount;
    private final Result mResult;

    /**
     * @param refCount the number of jobs that will make use of the result.
     * @param result the result of running a job.
     */
    CountResult(int refCount, Result result) {
      assert refCount > 0;
      mRefCount = refCount;
      mResult = result;
    }

    int decrement() {
      --mRefCount;
      assert mRefCount >= 0;
      return mRefCount;
    }

    int refCount() {
      return mRefCount;
    }

    Result result() {
      return mResult;
    }

    @Override
    public boolean integrity() {
      Exam.assertTrue(mRefCount >= 0);
      return true;
    }
  }

  private final Dependencies<J> mDependencies;

  private final JobFactory<J> mFactory;

  private final EventList<J> mEventList;

  private final PrintStream mTrace;

  private final Map<J, CountResult> mResults = new HashMap<>();

  private final Set<J> mRunning = new HashSet<>();

  private final JobStatistics<J> mStatistics;

  private final LookAhead mLookAhead;

  @Override
  public LookAhead lookAhead() {
    return mLookAhead;
  }

  @Override
  public boolean checkEmpty() {
    Exam.assertEquals(null, mEventList.next(mLookAhead));
    Exam.assertEquals(mResults.toString(), 0, mResults.size());
    Exam.assertEquals(0, mRunning.size());
    return true;
  }

  @Override
  public synchronized boolean integrity() {
    Exam.assertNotNull(mDependencies);
    Exam.assertNotNull(mResults);
    Exam.assertNotNull(mRunning);
    Exam.assertNotNull(mEventList);
    Exam.assertNotNull(mLookAhead);
    return true;
  }

  @Override
  public synchronized boolean globalIntegrity() {
    integrity();
    //results running and candidates should be disjoint
    for (final Map.Entry<J, CountResult> e : mResults.entrySet()) {
      final J id = e.getKey();
      Exam.assertTrue(e.getValue().refCount() > 0);
      Exam.assertFalse(mRunning.contains(id));
      Exam.assertFalse(mEventList.contains(id));
    }
    for (final J id : mRunning) {
      Exam.assertFalse(mEventList.contains(id));
    }
    return true;
  }

  /**
   * @param dependencies computes the dependencies between job identifiers and converts them to
   * @param factory converts {@link JobId}s to runnable jobs.
   * @param eventList event list to use for scheduling jobs that have their parameters available (will be tailored to the particular job types).
   * @param trace if non-null then trace the job identifiers of each job as they are started and finished.
   * @param statistics record statistics about execution of jobs (if non-null).
   * @param lookAhead number of chunks ahead of the earliest time that chunks can be scheduled.
   */
  public SchedulerSynchronized(final Dependencies<J> dependencies, final JobFactory<J> factory, final EventList<J> eventList, final PrintStream trace, final JobStatistics<J> statistics, final int lookAhead) {
    mDependencies = dependencies;
    mFactory = factory;
    mEventList = eventList;
    mTrace = trace;
    mStatistics = statistics;
    mLookAhead = new LookAhead(lookAhead, mDependencies.delta());
  }

  private final double[] mStarvationCounts = new double[JobType.values().length];

  @Override
  public synchronized Job<J> doneNext(J id, Result result, long nanoTime) {
    if (mStatistics != null && id != null) {
      mStatistics.increment(id, nanoTime);
    }
    trace(">", id, id == null ? "" : result);
    if (id != null) {
      assert id.validResult(result); // : "validating result id=" + id + " result=" + result;
      mRunning.remove(id);
      final Collection<J> to = mDependencies.to(id);
      assert Util.checkOrder(id, to, +1);
      final int size = Util.nonNullSize(to);
      if (size > 0) {
        mResults.put(id, new CountResult(size, result));
      }
      for (final J idTo : to) {
        if (idTo != null && argumentsAvailable(idTo)) {
          trace("*", idTo, "");
          mLookAhead.increment(idTo.time());
          mEventList.add(idTo);
        }
      }
      mLookAhead.decrement(id.time());
    }
    final Job<J> runnableJob = getRunnableJob();
    if (runnableJob != null) {
      mRunning.add(runnableJob.id());
    } else {
      //System.out.println("Starving: " + mRunning);
      final double inc = 1.0 / mRunning.size();
      for (final J job : mRunning) {
        if (job instanceof JobIdMultisample) {
          mStarvationCounts[((JobIdMultisample) job).type().ordinal()] += inc;
        }
      }
    }
    trace("<", runnableJob == null ? null : runnableJob.id(), "");
    return runnableJob;
  }

  /** Log statistics relating to starvation. */
  public void dumpStarvation() {
    if (ArrayUtils.sum(mStarvationCounts) > 0) {
      final double[] normalized = MathUtils.renormalize(mStarvationCounts);
      Diagnostic.developerLog("Starvation statisitcs:");
      for (final JobType jt : JobType.values()) {
        Diagnostic.developerLog(jt + " " + Utils.realFormat(100.0 * normalized[jt.ordinal()], 2) + "%");
      }
    }
  }

  /**
   * Trace execution of a job.
   * @param inOut indicates if being started or finished.
   * @param id identifier of job.
   * @param thing additional message.
   */
  private void trace(final String inOut, final J id, final Object thing) {
    if (mTrace == null) {
      return;
    }
    final String msg;
    if (thing == null) {
      msg = "null";
    } else if (thing instanceof String) {
      msg = (String) thing;
    } else if (thing.getClass().isArray()) {
      msg = Arrays.toString((Object[]) thing);
    } else {
      msg = thing.toString();
    }
    mTrace.print(inOut);
    if (id != null) {
      mTrace.print(id);
    }
    mTrace.println(" " + msg);
    mTrace.flush();
  }

  /**
   * check if id has all its arguments computed and so is a candidate to be executed.
   * @param id the job identifier being checked.
   * @return true if id may be executed
   */
  boolean argumentsAvailable(final J id) {
    final Collection<J> from = mDependencies.from(id);
    assert Util.checkOrder(id, from, -1);
    for (final J frId : from) {
      if (frId != null && !mResults.containsKey(frId)) {
        return false;
      }
    }
    return true;
  }

  private Job<J> getRunnableJob() {
    final J next0 = mEventList.next(mLookAhead);
    if (next0 != null) {
      return j2Job(next0);
    }
    final J next = mDependencies.next(mLookAhead);
    if (next == null) {
      return null;
    }
    mLookAhead.increment(next.time());
    return j2Job(next);
  }

  /**
   * Get a job to be run from <code>mDependencies</code> using the arguments stored in <code>mResults</code>.
   * Garbage collect <code>mResults</code> as all uses of an entry are satisfied.
   * @param id the job identifier to be converted.
   * @return the job.
   */
  Job<J> j2Job(final J id) {
    //get the from arguments
    final Collection<J> from = mDependencies.from(id);
    assert Util.checkOrder(id, from, -1);
    final Result[] arguments = new Result[from.size()];
    final Iterator<J> it = from.iterator();
    for (int i = 0; it.hasNext(); ++i) {
      final J argid = it.next();
      final Result argv;
      if (argid == null) {
        argv = null;
      } else {
        final CountResult countResult = mResults.get(argid);
        if (countResult.decrement() == 0) {
          mResults.remove(argid);
        }
        argv = countResult.result();
      }
      arguments[i] = argv;
    }
    final Job<J> job = mFactory.job(id, arguments);
    trace("+", id, arguments);
    assert id.validArguments(arguments); // : "validating arguments" + id + Arrays.toString(arguments);
    return job;
  }
}
