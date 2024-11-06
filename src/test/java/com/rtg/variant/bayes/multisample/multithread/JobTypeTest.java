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
import java.util.Collections;

import com.rtg.scheduler.Result;
import com.rtg.util.TestUtils;
import com.rtg.variant.Variant;
import com.rtg.variant.bayes.multisample.Complexities;

import junit.framework.TestCase;

/**
 */
public class JobTypeTest extends TestCase {

  public void test() {
    TestUtils.testEnum(JobType.class, "[INCR, DANGLING, COMPLEX, FLUSH, FILTER, BED, OUT]");
  }

  public void testIncrement() {
    final JobType t = JobType.INCR;
    assertTrue(t.validResult(getValidIncrementResult()));
    assertFalse(t.validArguments(new Result[] {new Result(new Object())}));
    assertTrue(t.validArguments(new Result[] {}));


    final Complexities regions = getComplexities();
    assertFalse(t.validResult(new Result(5, regions)));
  }

  public void testFlush() {
    final JobType t = JobType.FLUSH;
    final Result result = new Result();
    assertTrue(t.validResult(result));
    assertFalse(t.validArguments(new Result[] {}));
    assertFalse(t.validArguments(new Result[] {result}));
    assertFalse(t.validArguments(new Result[] {result, result, result}));
    assertFalse(t.validArguments(new Result[] {getValidIncrementResult(), result, result}));
    assertFalse(t.validArguments(new Result[] {result, getValidIncrementResult(), result}));
    assertFalse(t.validArguments(new Result[] {getValidIncrementResult(), getValidIncrementResult(), result, result}));
    assertTrue(t.validArguments(new Result[] {getValidIncrementResult(), getValidComplexResult()}));
  }

  private Complexities getComplexities() {
    final ArrayList<Variant> chunk = new ArrayList<>();
    chunk.add(com.rtg.variant.bayes.multisample.TestUtils.createVariant(4));
    chunk.add(com.rtg.variant.bayes.multisample.TestUtils.createVariant(5));
    chunk.add(com.rtg.variant.bayes.multisample.TestUtils.createVariant(6));
    chunk.add(com.rtg.variant.bayes.multisample.TestUtils.createVariant(10));
    chunk.add(com.rtg.variant.bayes.multisample.TestUtils.createVariant(15));

    return new Complexities(chunk, "foo", 0, 50, 3, 15, new byte[] {0, 1, 2, 3, 4, 0}, true, null);
  }

  public void testDangling() {
    final JobType t = JobType.DANGLING;
    assertTrue(t.validResult(getValidDanglingResult()));

    assertFalse(t.validArguments(new Result[] {new Result(new Object())}));
    final Result r = getValidIncrementResult();
    final Complexities regions2 = getComplexities();
    final Result r2 = new Result(regions2, 8);
    final Result r3 = getValidDanglingResult();
    assertTrue(t.validArguments(new Result[] {r, r2, r3}));
    assertTrue(t.validArguments(new Result[] {null, r2, null}));
    assertTrue(t.validArguments(new Result[] {r, null, r3}));
    assertTrue(t.validArguments(new Result[] {r, r2, null}));
    assertFalse(t.validArguments(new Result[] {r, new Result(), r3}));
    assertFalse(t.validArguments(new Result[] {new Result(), r, r3}));
    assertFalse(t.validResult(new Result()));
    assertTrue(t.validResult(new Result((Object) null)));
  }

  public void testComplex() {
    final JobType t = JobType.COMPLEX;
    assertTrue(t.validResult(getValidComplexResult()));
    assertFalse(t.validArguments(new Result[] {new Result(new Object())}));
    assertTrue(t.validArguments(new Result[] {new Result((Object) null)}));

    final Result r = getValidDanglingResult();
    assertTrue(t.validArguments(new Result[] {r}));
    assertFalse(t.validArguments(new Result[] {new Result()}));
    assertTrue(t.validResult(new Result(null, null)));
    assertFalse(t.validResult(new Result()));
  }

  public void testFilter() {
    final JobType t = JobType.FILTER;
    assertTrue(t.validResult(getValidFilterResult()));
    assertFalse(t.validArguments(new Result[] {new Result(new Object())}));

    final ArrayList<Variant> list = new ArrayList<>();
    assertTrue(t.validArguments(new Result[] {null, null}));

    final Result r = getValidIncrementResult(); //increment
    final Result r2 = getValidComplexResult(); //complex
    assertTrue(t.validArguments(new Result[] {r, r2}));
    assertFalse(t.validArguments(new Result[] {new Result(), r2}));
    assertFalse(t.validArguments(new Result[] {r, new Result()}));

    assertFalse(t.validResult(new Result()));
    assertTrue(t.validResult(new Result(null, 42)));
    assertTrue(t.validResult(new Result(list, null)));
    assertTrue(t.validResult(new Result(null, null)));
  }

  public void testBed() {
    final JobType t = JobType.BED;
    assertTrue(t.validResult(getValidBedResult()));
    assertFalse(t.validArguments(new Result[] {new Result(new Object())}));
    final Result r1 = getValidBedResult();
    final Result r3 = getValidComplexResult();
    assertTrue(t.validArguments(new Result[] {null, r3})); //for time 0
    assertTrue(t.validArguments(new Result[] {r1, r3}));
    assertFalse(t.validArguments(new Result[] {r3, r1}));
    assertFalse(t.validResult(new Result(new Object())));
  }

  public void testOut() {
    final JobType t = JobType.OUT;
    assertTrue(t.validResult(getValidOutResult()));
    assertFalse(t.validArguments(new Result[] {new Result(new Object())}));

    final Result r1 = getValidOutResult();
    final Result r2 = getValidFilterResult();

    assertTrue(t.validArguments(new Result[] {r1, r2, new Result()}));
    assertFalse(t.validResult(new Result()));
  }

  private Result getValidIncrementResult() {
    final Complexities regions = getComplexities();
    return new Result(regions, 5);
  }

  private Result getValidComplexResult() {
    final ArrayList<Variant> list = new ArrayList<>();
    return new Result(list, getComplexities());
  }

  private Result getValidFilterResult() {
    final ArrayList<Variant> list = new ArrayList<>();
    return new Result(list, 42);
  }

  private Result getValidBedResult() {
    return new Result();
  }

  private Result getValidDanglingResult() {
    final Complexities regions = getComplexities();
    return new Result(regions);
  }

  private Result getValidOutResult() {
    return new Result(Collections.emptyList());
  }

}

