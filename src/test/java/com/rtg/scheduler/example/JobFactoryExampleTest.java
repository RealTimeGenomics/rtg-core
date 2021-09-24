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

import com.rtg.scheduler.Job;
import com.rtg.scheduler.JobFactory;

import junit.framework.TestCase;


/**
 */
public class JobFactoryExampleTest extends TestCase {

  private ExampleJobId id(final int ch, final JobType ty) {
    return new ExampleJobId(5, ch, ty);
  }

  public void testJob() {
    final JobFactory<ExampleJobId> jf = new JobFactoryExample(5);
    final Job<ExampleJobId> job = jf.job(id(2, JobType.INCR), new MockResult[] {new MockResult("arg0"), new MockResult("arg1")});
    assertEquals("2:INCR(arg0, arg1)", job.toString());
  }
}
