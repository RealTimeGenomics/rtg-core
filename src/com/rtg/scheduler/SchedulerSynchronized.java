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

package com.rtg.scheduler;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.util.integrity.Integrity;

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
    public CountResult(int refCount, Result result) {
      assert refCount > 0;
      mRefCount = refCount;
      mResult = result;
    }

    int decrement() {
      mRefCount--;
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
    for (final J id : mResults.keySet()) {
      Exam.assertTrue(mResults.get(id).refCount() > 0);
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
    }
    trace("<", runnableJob == null ? null : runnableJob.id(), "");
    return runnableJob;
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
    for (int i = 0; it.hasNext(); i++) {
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
