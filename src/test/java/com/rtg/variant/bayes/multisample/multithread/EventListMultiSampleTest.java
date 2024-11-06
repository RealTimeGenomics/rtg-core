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

package com.rtg.variant.bayes.multisample.multithread;

import com.rtg.scheduler.LookAhead;
import com.rtg.scheduler.enumtime.EnumTimeId;

import junit.framework.TestCase;

/**
 */
public class EventListMultiSampleTest extends TestCase {

  private static class MockLookAhead extends LookAhead {
    final boolean mOk;

    MockLookAhead(boolean ok) {
      super(1, 0);
      mOk = ok;
    }
    @Override
    public synchronized boolean ok(int t, int delta) {
      assertEquals(0, delta);
      return mOk;
    }
  }

  public void test() {
    final EventListMultiSample<JobIdMultisample> ev = new EventListMultiSample<>();
    ev.globalIntegrity();
    final JobIdMultisample id1 = new JobIdMultisample(5, 0, JobType.INCR);
    assertFalse(ev.contains(id1));
    ev.add(id1);
    assertTrue(ev.contains(id1));
    //System.err.println(ev);
    ev.globalIntegrity();
    final EnumTimeId<?> n1 = ev.next(new MockLookAhead(true));
    assertTrue(n1.time() == 0 && n1.type() == JobType.INCR);
    ev.globalIntegrity();
    assertNull(ev.next(new MockLookAhead(true)));

    final JobIdMultisample id2 = new JobIdMultisample(5, 2, JobType.INCR);
    assertFalse(ev.contains(id2));
    ev.add(id2);
    assertTrue(ev.contains(id2));
    ev.add(new JobIdMultisample(5, 2, JobType.FILTER));
    ev.add(new JobIdMultisample(5, 0, JobType.INCR));
    //System.err.println(ev);
    ev.globalIntegrity();
    final EnumTimeId<?> n2 = ev.next(new MockLookAhead(true));
    assertTrue(n2.time() == 0 && n2.type() == JobType.INCR);
    final EnumTimeId<?> n3 = ev.next(new MockLookAhead(true));
    assertTrue(n3.time() == 2 && n3.type() == JobType.INCR);
    assertNull(ev.next(new MockLookAhead(false)));
    final EnumTimeId<?> n4 = ev.next(new MockLookAhead(true));
    assertTrue(n4.time() == 2 && n4.type() == JobType.FILTER);
    assertNull(ev.next(new MockLookAhead(true)));

    ev.add(new JobIdMultisample(5, 1, JobType.INCR));
    ev.add(new JobIdMultisample(5, 2, JobType.FILTER));
    ev.add(new JobIdMultisample(5, 0, JobType.INCR));
    ev.globalIntegrity();
  }

  public void testWrap() {
    final EventListMultiSample<JobIdMultisample> ev = new EventListMultiSample<>();
    ev.globalIntegrity();
    ev.add(new JobIdMultisample(5, 2, JobType.INCR));
    ev.add(new JobIdMultisample(5, 3, JobType.INCR));
    ev.add(new JobIdMultisample(5, 4, JobType.INCR));
    ev.globalIntegrity();
    final EnumTimeId<?> n1 = ev.next(new MockLookAhead(true));
    assertTrue(n1.time() == 2 && n1.type() == JobType.INCR);
    ev.globalIntegrity();
    final EnumTimeId<?> n2 = ev.next(new MockLookAhead(true));
    assertTrue(n2.time() == 3 && n2.type() == JobType.INCR);
    ev.globalIntegrity();
    final EnumTimeId<?> n3 = ev.next(new MockLookAhead(true));
    assertTrue(n3.time() == 4 && n3.type() == JobType.INCR);
    ev.globalIntegrity();
    assertNull(ev.next(new MockLookAhead(true)));
  }

  public void testNewLength() {
    assertEquals(3, EventListMultiSample.newLength(1, 2));
    assertEquals(2, EventListMultiSample.newLength(1, 1));
    assertEquals(12, EventListMultiSample.newLength(3, 10));
  }
}
