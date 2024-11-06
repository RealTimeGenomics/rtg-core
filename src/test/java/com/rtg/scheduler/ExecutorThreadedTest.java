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

import java.io.IOException;
import java.util.concurrent.atomic.AtomicInteger;

import com.rtg.scheduler.example.DependenciesExample;
import com.rtg.scheduler.example.ExampleJobId;
import com.rtg.scheduler.example.JobFactoryExample;
import com.rtg.scheduler.example.JobType;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.variant.bayes.multisample.multithread.EventListMultiSample;

import junit.framework.TestCase;

/**
 */
public class ExecutorThreadedTest extends TestCase {

  //small test with explicit answers
  public void test2() throws IOException {
    Diagnostic.setLogStream();
    final DependenciesExample dep = new DependenciesExample(2);
    final JobFactoryExample factory = new JobFactoryExample(2);
    final MemoryPrintStream ps = new MemoryPrintStream();
    final EventList<ExampleJobId> eventList = new EventListMultiSample<>();
    final Scheduler<ExampleJobId> sch = new SchedulerSynchronized<>(dep, factory, eventList, ps.printStream(), null, 2);
    new ExecutorThreaded<>(sch, 5).run();
    assertNull(eventList.next(sch.lookAhead()));
    final String str = ps.toString();
    //System.err.println(str);
    TestUtils.containsAll(str, ">", "<", "<0:INCR", ">0:INCR", "<1:OUT", ">1:OUT");
    assertEquals("1:BED(0:BED(null, 0:DANGLING(null, 0:INCR())), 1:DANGLING(0:INCR(), 1:INCR()))", factory.bed().toString());
    assertEquals("1:OUT(0:OUT(null, 0:MERGE(0:DANGLING(null, 0:INCR()), 0:COMPLEX(0:DANGLING(null, 0:INCR())))), 1:MERGE(1:DANGLING(0:INCR(), 1:INCR()), 1:COMPLEX(1:DANGLING(0:INCR(), 1:INCR()))))",
        factory.out().toString());
  }

  //a larger example check that the results for sequential and random are the same
  public void test5() throws IOException {
    Diagnostic.setLogStream();
    final DependenciesExample depRand = new DependenciesExample(5);
    final JobFactoryExample factoryRand = new JobFactoryExample(5);
    final MemoryPrintStream psRand = new MemoryPrintStream();
    final EventList<ExampleJobId> eventList = new EventListMultiSample<>();
    final Scheduler<ExampleJobId> schRand = new SchedulerSynchronized<>(depRand, factoryRand, eventList, psRand.printStream(), null, 5);
    new ExecutorThreaded<>(schRand, 5).run();
    assertNull(eventList.next(schRand.lookAhead()));
    final String strRand = psRand.toString();
    //System.err.println(strRand);
    TestUtils.containsAll(strRand, ">", "<", "<0:INCR", ">0:INCR", "<4:OUT", ">4:OUT");
    final String bedRand = factoryRand.bed().toString();
    final String outRand = factoryRand.out().toString();

    final DependenciesExample depSeq = new DependenciesExample(5);
    final JobFactoryExample factorySeq = new JobFactoryExample(5);
    final MemoryPrintStream psSeq = new MemoryPrintStream();
    final Scheduler<ExampleJobId> schSeq = new SchedulerSynchronized<>(depSeq, factorySeq, eventList, psSeq.printStream(), null, 5);
    new ExecutorSequential<>(schSeq).run();
    assertNull(eventList.next(schSeq.lookAhead()));
    final String strSeq = psSeq.toString();
    //System.err.println(strSeq);
    TestUtils.containsAll(strSeq, ">", "<", "<0:INCR", ">0:INCR", "<4:OUT", ">4:OUT");
    final String bedSeq = factorySeq.bed().toString();
    final String outSeq = factorySeq.out().toString();

    assertFalse(strRand.equals(strSeq));
    assertEquals(bedRand, bedSeq);
    assertEquals(outRand, outSeq);
  }

  class MockJobError extends Job<ExampleJobId> {

    MockJobError(ExampleJobId id) {
      super(id);
    }

    @Override
    protected Result run() {
      throw new RuntimeException("KILL");
    }

  }

  class MockJobNoError extends Job<ExampleJobId> {

    private final AtomicInteger mAtom;

    MockJobNoError(ExampleJobId id, final AtomicInteger atomI) {
      super(id);
      mAtom = atomI;
    }

    @Override
    protected Result run() {
      mAtom.incrementAndGet();
      return new Result();
    }
  }

  class MockSched implements Scheduler<ExampleJobId> {

    private int mStart = 0;
    private final AtomicInteger mAtom;

    MockSched(final AtomicInteger i) {
      mAtom = i;
    }

    @Override
    public synchronized Job<ExampleJobId> doneNext(ExampleJobId id, Result result, long nanoTime) {
      if (mStart > 9) {
        return null;
      }
      // when this reaches 3, fire the job so we get an error
      // after which no other Jobs should run
      if (mAtom.get() >= 3) {
        return new MockJobError(new ExampleJobId(10, mStart++, JobType.INCR));
      }
      return new MockJobNoError(new ExampleJobId(10, mStart++, JobType.INCR), mAtom);
    }

    @Override
    public LookAhead lookAhead() {
      throw new UnsupportedOperationException();
    }

    @Override
    public boolean checkEmpty() {
      return true;
    }
  }

  public void testEarlyExit() {
    Diagnostic.setLogStream();
    final AtomicInteger atom = new AtomicInteger();
    atom.set(0);
    final ExecutorThreaded<ExampleJobId> t = new ExecutorThreaded<>(new MockSched(atom), 2);
    try {
      t.run();
      fail();
    } catch (final IOException e) {
    } catch (final RuntimeException e) {
      assertEquals("KILL", e.getCause().getMessage());
    }
  }
}
