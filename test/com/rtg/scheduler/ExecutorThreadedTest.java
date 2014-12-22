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

    public MockJobError(ExampleJobId id) {
      super(id);
    }

    @Override
    protected Result run() {
      throw new RuntimeException("KILL");
    }

  }

  class MockJobNoError extends Job<ExampleJobId> {

    private final AtomicInteger mAtom;

    public MockJobNoError(ExampleJobId id, final AtomicInteger atomI) {
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

    public MockSched(final AtomicInteger i) {
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
