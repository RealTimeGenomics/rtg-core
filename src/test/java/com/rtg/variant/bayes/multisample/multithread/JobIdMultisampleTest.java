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
