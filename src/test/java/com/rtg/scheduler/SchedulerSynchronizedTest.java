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
