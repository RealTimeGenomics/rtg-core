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

package com.rtg.scheduler.enumtime;

import java.util.Collection;

import com.rtg.scheduler.example.DependenciesExample;
import com.rtg.scheduler.example.ExampleJobId;
import com.rtg.scheduler.example.JobType;
import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class LocalOrderingTest extends TestCase {

  private ExampleJobId id(final int ch, final JobType ty) {
    return new ExampleJobId(5, ch, ty);
  }

  public void testTo() {
    final DependenciesExample ordering = new DependenciesExample(5);
    checkTo(ordering, JobType.INCR, id(0, JobType.DANGLING), id(1, JobType.DANGLING));
    checkTo(ordering, JobType.DANGLING, id(0, JobType.BED), id(0, JobType.COMPLEX), id(0, JobType.MERGE));
    checkTo(ordering, JobType.BED, id(1, JobType.BED));
    checkTo(ordering, JobType.COMPLEX, id(0, JobType.MERGE));
    checkTo(ordering, JobType.MERGE, id(0, JobType.OUT));
    checkTo(ordering, JobType.OUT, id(1, JobType.OUT));
  }

  private void checkTo(final DependenciesExample dep, final JobType ty, final ExampleJobId... exp) {
    final ExampleJobId id = id(0, ty);
    final Collection<ExampleJobId> to = dep.to(id);
    assertEquals(exp.length, to.size());
    for (final ExampleJobId exid : exp) {
      to.contains(exid);
      assertTrue(id + ">" + to.toString(), to.contains(exid));
    }
  }

  public void testFrom() {
    final DependenciesExample dep = new DependenciesExample(5);
    checkFrom(dep, JobType.INCR);
    checkFrom(dep, JobType.DANGLING, id(0, JobType.INCR), id(1, JobType.INCR));
    checkFrom(dep, JobType.BED, id(0, JobType.BED), id(1, JobType.DANGLING));
    checkFrom(dep, JobType.COMPLEX, id(1, JobType.DANGLING));
    checkFrom(dep, JobType.MERGE, id(1, JobType.DANGLING), id(1, JobType.COMPLEX));
    checkFrom(dep, JobType.OUT, id(0, JobType.OUT), id(1, JobType.MERGE));
  }

  private void checkFrom(final DependenciesExample dep, final JobType ty, final ExampleJobId... exp) {
    final ExampleJobId id = id(1, ty);
    final Collection<ExampleJobId> from = dep.from(id);
    assertEquals(exp.length, from.size());
    int i = 0;
    for (final Object frid : from) {
      assertTrue(id + ">" + from.toString(), exp[i].equals(frid));
      ++i;
    }
  }

  public void testOrdering() {
    final ExampleJobId id0 = new ExampleJobId(5, 0, JobType.INCR);
    final ExampleJobId id1 = new ExampleJobId(5, 1, JobType.INCR);
    final ExampleJobId id2 = new ExampleJobId(5, 1, JobType.DANGLING);
    final ExampleJobId id3a = new ExampleJobId(5, 1, JobType.COMPLEX);
    final ExampleJobId id3b = new ExampleJobId(5, 1, JobType.BED);
    final ExampleJobId id4 = new ExampleJobId(5, 1, JobType.MERGE);
    final ExampleJobId id5 = new ExampleJobId(5, 1, JobType.OUT);
    TestUtils.testOrder(new ExampleJobId[] {id0, id1, id2, id3a, id4, id5}, true);
    TestUtils.testOrder(new ExampleJobId[] {id0, id1, id2, id3b}, true);
  }

  public void testFrozen() {
    final LocalOrdering<JobType> ord = new LocalOrdering<>(JobType.values());
    assertFalse(ord.frozen());
    try {
      ord.to(JobType.INCR);
      fail();
    } catch (final RuntimeException e) {
      //Expected
    }
    try {
      ord.from(JobType.INCR);
      fail();
    } catch (final RuntimeException e) {
      //Expected
    }
    try {
      ord.before(JobType.INCR, JobType.INCR);
      fail();
    } catch (final RuntimeException e) {
      //Expected
    }
    ord.freeze();
    assertTrue(ord.frozen());
    try {
      ord.setLink(JobType.INCR, 0, JobType.COMPLEX);
      fail();
    } catch (final RuntimeException e) {
      //Expected
    }
  }

  public void testTransitiveClosure1() {
    final boolean[][] m = new boolean[3][3];
    m[1][0] = true;
    m[0][2] = true;
    LocalOrdering.transitiveClosure(m);
    assertTrue(m[1][2]);
  }

  public void testTransitiveClosure2() {
    final boolean[][] m = new boolean[3][3];
    m[2][1] = true;
    m[1][0] = true;
    LocalOrdering.transitiveClosure(m);
    assertTrue(m[2][0]);
  }
}
