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
package com.rtg.util;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;


/**
 * A daemon thread that can be given more than one <code>Runnable</code> over
 * its lifetime. Used by <code>SimpleThreadPool</code>
 */
public class WorkerThread extends Thread {

  private static final int SLEEP_TIME = 5 * 1000; //5 seconds
  private final Object mSleepSync = new Object();
  private final Object mJobSync = new Object();

  private final Object mCompleteNotify;

  private boolean mDie = false;

  private volatile Runnable mJob;

  /**
   * Constructs a worker thread
   * @param name name of thread
   * @param complete Object on which to call <code>notifyAll</code> on when a task is completed
   */
  public WorkerThread(final String name, final Object complete) {
    super(name);
    setDaemon(true);
    mCompleteNotify = complete;
  }

  /**
   * Enqueue a job for this thread
   * @param job ob to run
   */
  public void enqueueJob(final Runnable job) {
    synchronized (mSleepSync) { //prevent race condition
      if (hasJob()) {
        throw new IllegalStateException("Job already enqueued");
      }
      setJob(job);
      //interrupt();
      mSleepSync.notifyAll();
    }
  }

  /**
   * Kill the thread
   */
  public void die() {
    synchronized (mSleepSync) {
      mDie = true;
      mSleepSync.notifyAll();
    }
  }

  private void setJob(final Runnable job) {
    synchronized (mJobSync) {
      mJob = job;
    }
  }

  /**
   * Checks whether the job queue is empty
   * @return true if empty
   */
  public boolean hasJob() {
    return mJob != null;
  }

  /**
   * Runs jobs in the job queue
   */
  @Override
  public void run() {
    while (!mDie) {
      synchronized (mJobSync) {
        if (hasJob()) {
          final Runnable job = mJob;
          try {
            job.run();
          } catch (final Throwable t) {
            Diagnostic.error(ErrorType.SLIM_ERROR);
            Diagnostic.userLog(t);
          }
          synchronized (mCompleteNotify) {
            setJob(null);
            mCompleteNotify.notifyAll();
          }
        }
      }
      synchronized (mSleepSync) {
        try {
          if (!hasJob() && !mDie) {
            mSleepSync.wait(SLEEP_TIME);
          }
        } catch (final InterruptedException e) {
          //no problem
        }
      }
    }
  }
}
