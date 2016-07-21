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
    assertTrue(t.validArguments(new Result[] {getValidIncrementResult(), getValidIncrementResult(), result}));
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

    assertTrue(t.validResult(new Result((Object) null)));
    assertFalse(t.validResult(new Result()));
  }

  public void testFilter() {
    final JobType t = JobType.FILTER;
    assertTrue(t.validResult(getValidFilterResult()));
    assertFalse(t.validArguments(new Result[] {new Result(new Object())}));

    final ArrayList<Variant> list = new ArrayList<>();
    assertTrue(t.validArguments(new Result[] {null, null, new Result(null, list), null}));
    assertTrue(t.validArguments(new Result[] {null, null, new Result(list, null), null}));


    final Result r = getValidIncrementResult(); //increment
    final Result r2 = getValidComplexResult(); //complex
    final Result r3 = getValidFilterResult(); //filter
    assertTrue(t.validArguments(new Result[] {r, r2, r3, new Result()}));
    assertTrue(t.validArguments(new Result[] {r, r2, null, new Result()}));
    assertFalse(t.validArguments(new Result[] {new Result(), r2, null, new Result()}));
    assertFalse(t.validArguments(new Result[] {r, new Result(), null, new Result()}));

    assertFalse(t.validResult(new Result()));
    assertTrue(t.validResult(new Result(null, list)));
    assertTrue(t.validResult(new Result(list, null)));
    assertTrue(t.validResult(new Result(null, null)));
  }

  public void testBed() {
    final JobType t = JobType.BED;
    assertTrue(t.validResult(getValidBedResult()));
    assertFalse(t.validArguments(new Result[] {new Result(new Object())}));
    final Result r1 = getValidBedResult();
    final Result r2 = getValidDanglingResult();
    final Result r3 = getValidComplexResult();
    assertTrue(t.validArguments(new Result[] {null, r2, r3})); //for time 0
    assertTrue(t.validArguments(new Result[] {r1, r2, r3}));
    assertFalse(t.validArguments(new Result[] {r2, r3, r1}));
    assertFalse(t.validResult(new Result(new Object())));
  }

  public void testOut() {
    final JobType t = JobType.OUT;
    assertTrue(t.validResult(getValidOutResult()));
    assertFalse(t.validArguments(new Result[] {new Result(new Object())}));

    final Result r1 = getValidOutResult();
    final Result r2 = getValidFilterResult();

    assertTrue(t.validArguments(new Result[] {r1, r2}));

    assertFalse(t.validResult(new Result()));
  }

  private Result getValidIncrementResult() {
    final Complexities regions = getComplexities();
    return new Result(regions, 5);
  }

  private Result getValidComplexResult() {
    final ArrayList<Variant> list = new ArrayList<>();
    return new Result(list);
  }

  private Result getValidFilterResult() {
    final ArrayList<Variant> list = new ArrayList<>();
    return new Result(list, list);
  }

  private Result getValidBedResult() {
    return new Result();
  }

  private Result getValidDanglingResult() {
    final Complexities regions = getComplexities();
    return new Result(regions);
  }

  private Result getValidOutResult() {
    return null;
  }

}

