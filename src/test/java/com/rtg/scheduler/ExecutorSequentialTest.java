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

import com.rtg.scheduler.example.DependenciesExample;
import com.rtg.scheduler.example.ExampleJobId;
import com.rtg.scheduler.example.JobFactoryExample;
import com.rtg.util.TestUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.variant.bayes.multisample.multithread.EventListMultiSample;

import junit.framework.TestCase;

/**
 */
public class ExecutorSequentialTest extends TestCase {

  public void test() throws IOException {
    final DependenciesExample dep = new DependenciesExample(2);
    final JobFactoryExample factory = new JobFactoryExample(2);
    final MemoryPrintStream ps = new MemoryPrintStream();
    final EventList<ExampleJobId> eventList = new EventListMultiSample<>();
    final Scheduler<ExampleJobId> sch = new SchedulerSynchronized<>(dep, factory, eventList, ps.printStream(), null, 2);
    new ExecutorSequential<>(sch).run();
    assertNull(eventList.next(sch.lookAhead()));
    final String str = ps.toString();
    TestUtils.containsAll(str, ">", "<", "<0:INCR", ">0:INCR", "<1:OUT", ">1:OUT");
    assertEquals("1:BED(0:BED(null, 0:DANGLING(null, 0:INCR())), 1:DANGLING(0:INCR(), 1:INCR()))", factory.bed().toString());
    assertEquals("1:OUT(0:OUT(null, 0:MERGE(0:DANGLING(null, 0:INCR()), 0:COMPLEX(0:DANGLING(null, 0:INCR())))), 1:MERGE(1:DANGLING(0:INCR(), 1:INCR()), 1:COMPLEX(1:DANGLING(0:INCR(), 1:INCR()))))",
        factory.out().toString());
  }
}
