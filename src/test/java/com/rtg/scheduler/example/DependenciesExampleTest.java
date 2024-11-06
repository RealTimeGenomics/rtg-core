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

package com.rtg.scheduler.example;

import java.util.Collection;
import java.util.Set;

import com.rtg.scheduler.AbstractDependenciesTest;
import com.rtg.scheduler.Dependencies;
import com.rtg.util.integrity.Exam;

/**
 */
public class DependenciesExampleTest extends AbstractDependenciesTest<ExampleJobId> {

  @Override
  protected Dependencies<ExampleJobId> getDependencies() {
    return new DependenciesExample(5);
  }

  public void testGlobalIntegrity() {
    Exam.globalIntegrity(getDependencies());
  }

  public void testDetails() {
    final Dependencies<ExampleJobId> dep = getDependencies();
    final Set<ExampleJobId> allIds = allIds(dep, 5);
    //System.err.println(allIds);
    final int numberChunks = 5;
    final JobType[] values = JobType.values();
    for (int i = 0; i < numberChunks; ++i) {
      for (final JobType ty : values) {
        assertTrue("i=" + i + " ty=" + ty, allIds.contains(new ExampleJobId(numberChunks, i, ty)));
      }
    }
  }

  private ExampleJobId id(final int ch, final JobType ty) {
    return new ExampleJobId(5, ch, ty);
  }

  public void testTo() {
    final Dependencies<ExampleJobId> dep = getDependencies();
    checkTo(dep, JobType.INCR, id(0, JobType.DANGLING), id(1, JobType.DANGLING));
    checkTo(dep, JobType.DANGLING, id(0, JobType.BED), id(0, JobType.COMPLEX), id(0, JobType.MERGE));
    checkTo(dep, JobType.BED, id(1, JobType.BED));
    checkTo(dep, JobType.COMPLEX, id(0, JobType.MERGE));
    checkTo(dep, JobType.MERGE, id(0, JobType.OUT));
    checkTo(dep, JobType.OUT, id(1, JobType.OUT));
  }

  private void checkTo(final Dependencies<ExampleJobId> dep, final JobType ty, final ExampleJobId... exp) {
    final ExampleJobId id = id(0, ty);
    final Collection<ExampleJobId> to = dep.to(id);
    assertEquals(exp.length, to.size());
    for (final ExampleJobId exid : exp) {
      assertTrue(id + ">" + to.toString(), to.contains(exid));
    }
  }

  public void testFrom() {
    final Dependencies<ExampleJobId> dep = getDependencies();
    checkFrom(dep, JobType.INCR);
    checkFrom(dep, JobType.DANGLING, id(0, JobType.INCR), id(1, JobType.INCR));
    checkFrom(dep, JobType.BED, id(0, JobType.BED), id(1, JobType.DANGLING));
    checkFrom(dep, JobType.COMPLEX, id(1, JobType.DANGLING));
    checkFrom(dep, JobType.MERGE, id(1, JobType.DANGLING), id(1, JobType.COMPLEX));
    checkFrom(dep, JobType.OUT, id(0, JobType.OUT), id(1, JobType.MERGE));
  }

  private void checkFrom(final Dependencies<ExampleJobId> dep, final JobType ty, final ExampleJobId... exp) {
    final ExampleJobId id = id(1, ty);
    final Collection<ExampleJobId> from = dep.from(id);
    assertEquals(exp.length, from.size());
    int i = 0;
    for (final Object frid : from) {
      assertTrue(id + ">" + from.toString(), exp[i].equals(frid));
      ++i;
    }
  }
}
