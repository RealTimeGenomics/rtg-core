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

import java.util.ArrayList;

import com.rtg.scheduler.Result;
import com.rtg.variant.bayes.multisample.Complexities;
import com.rtg.variant.bayes.multisample.ComplexitiesTest;

import junit.framework.TestCase;

/**
 */
public class JobIdMultisampleTest extends TestCase {

  public void test() {
    final JobIdMultisample jbs = new JobIdMultisample(20, 5, JobType.INCR);
    final Complexities complexities = new Complexities(new ArrayList<>(), null, 1, 2, 0, 0, ComplexitiesTest.template(40), true, null);
    final Result r = new Result(complexities, 5);
    assertTrue(jbs.validResult(r));
    final Result[] res = {};
    assertTrue(jbs.validArguments(res));
  }

  public void testCompare() {
    final JobIdMultisample jbs = new JobIdMultisample(20, 5, JobType.INCR);
    final JobIdMultisample jbs2 = new JobIdMultisample(20, 5, JobType.INCR);

    assertEquals(jbs, jbs2);
    assertEquals(0, jbs.compareTo(jbs2));
    assertEquals(jbs.hashCode(), jbs2.hashCode());
  }

}
