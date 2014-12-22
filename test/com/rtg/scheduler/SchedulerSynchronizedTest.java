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

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;

import com.rtg.scheduler.example.DependenciesExample;
import com.rtg.scheduler.example.ExampleJobId;
import com.rtg.scheduler.example.JobFactoryExample;
import com.rtg.util.TestUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.variant.bayes.multisample.multithread.EventListMultiSample;

import junit.framework.TestCase;

/**
 */
public class SchedulerSynchronizedTest extends TestCase {

  private static class MockStatistics implements JobStatistics<ExampleJobId> {
    private final StringBuilder mSb = new StringBuilder();

    @Override
    public void increment(ExampleJobId id, long nanoTime) {
      if (id == null) {
        assertEquals(-1, nanoTime);
      } else {
        assertTrue(nanoTime >= 0);
      }
      mSb.append(id).append(" ").append(nanoTime).append(LS);
    }

    @Override
    public String toString() {
      return mSb.toString();
    }

  }

  //small test with explicit answers
  public void test2() throws IOException {
    final DependenciesExample dep = new DependenciesExample(2);
    final JobFactoryExample factory = new JobFactoryExample(2);
    final MemoryPrintStream ps = new MemoryPrintStream();
    final MockStatistics st = new MockStatistics();
    final EventList<ExampleJobId> eventList = new EventListMultiSample<>();
    final Scheduler<ExampleJobId> sch = new SchedulerSynchronized<>(dep, factory, eventList, ps.printStream(), st, 2);
    sch.checkEmpty();
    new ExecutorSequential<>(sch).run();
    sch.checkEmpty();
    final String str = ps.toString();
    //System.err.println(str);
    TestUtils.containsAll(str, "> " + LS, "<", "<0:INCR", ">0:INCR 0:INCR()", "<1:OUT", ">1:OUT null", "+0:DANGLING [null, 0:INCR()]", "*0:OUT");
    assertEquals("1:BED(0:BED(null, 0:DANGLING(null, 0:INCR())), 1:DANGLING(0:INCR(), 1:INCR()))", factory.bed().toString());
    assertEquals("1:OUT(0:OUT(null, 0:MERGE(0:DANGLING(null, 0:INCR()), 0:COMPLEX(0:DANGLING(null, 0:INCR())))), 1:MERGE(1:DANGLING(0:INCR(), 1:INCR()), 1:COMPLEX(1:DANGLING(0:INCR(), 1:INCR()))))",
        factory.out().toString());

    final String stStr = st.toString();
    //System.err.println(stStr);
    assertTrue(stStr.matches("(?m)(?s).*^0:INCR \\d+$.*")); //WOOHOO
  }

  //a larger example check that the results for sequential and random are the same
  public void test5() throws IOException {
    final DependenciesExample depRand = new DependenciesExample(5);
    final JobFactoryExample factoryRand = new JobFactoryExample(5);
    final MemoryPrintStream psRand = new MemoryPrintStream();
    final EventList<ExampleJobId> eventList = new EventListMultiSample<>();
    final Scheduler<ExampleJobId> schRand = new SchedulerSynchronized<>(depRand, factoryRand, eventList, psRand.printStream(), null, 5);
    schRand.checkEmpty();
    new ExecutorRandom<>(schRand, 5, (long) 42).run();
    schRand.checkEmpty();
    final String strRand = psRand.toString();
    //System.err.println(strRand);
    TestUtils.containsAll(strRand, ">", "<", "<0:INCR", ">0:INCR", "<4:OUT", ">4:OUT");
    final String bedRand = factoryRand.bed().toString();
    final String outRand = factoryRand.out().toString();

    final DependenciesExample depSeq = new DependenciesExample(5);
    final JobFactoryExample factorySeq = new JobFactoryExample(5);
    final MemoryPrintStream psSeq = new MemoryPrintStream();
    final Scheduler<ExampleJobId> schSeq = new SchedulerSynchronized<>(depSeq, factorySeq, eventList, psSeq.printStream(), null, 5);
    schSeq.checkEmpty();
    new ExecutorSequential<>(schSeq).run();
    schSeq.checkEmpty();
    final String strSeq = psSeq.toString();
    //System.err.println(strSeq);
    TestUtils.containsAll(strSeq, ">", "<", "<0:INCR", ">0:INCR", "<4:OUT", ">4:OUT");
    final String bedSeq = factorySeq.bed().toString();
    final String outSeq = factorySeq.out().toString();

    assertFalse(strRand.equals(strSeq));
    assertEquals(bedRand, bedSeq);
    assertEquals(outRand, outSeq);
  }

}
