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
    for (int i = 0; i < numberChunks; i++) {
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
      i++;
    }
  }
}
