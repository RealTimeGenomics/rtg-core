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

import java.util.concurrent.atomic.AtomicInteger;

import junit.framework.TestCase;

/**
 *
 */
public class WorkerThreadTest extends TestCase {
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(WorkerThreadTest.class);
  }

  public void testThread() {
    final AtomicInteger val = new AtomicInteger();
    val.set(0);
    final Runnable run = new Runnable() {
      @Override
      public void run() {
        val.incrementAndGet();
      }
    };
    final Object signal = new Object();
    final WorkerThread wt = new WorkerThread("testWorkerThread", signal);
    assertTrue(wt.isDaemon()); // sigh jumble
    wt.start();
    for (int i = 0; i < 10000; i++) {
      synchronized (signal) {
        wt.enqueueJob(run);
        while (wt.hasJob()) {
          try {
            signal.wait(1000);
          } catch (final InterruptedException e) {
            //
          }
        }
      }
    }
    wt.die();
    assertEquals(10000, val.get());
    try {
      Thread.sleep(100);
    } catch (final InterruptedException e) {

    }
    assertTrue(!wt.isAlive());
  }
}
