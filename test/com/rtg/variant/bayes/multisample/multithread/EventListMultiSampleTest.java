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
