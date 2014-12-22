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

package com.rtg.scheduler;

import com.rtg.util.integrity.IntegralAbstract;

import junit.framework.TestCase;

/**
 */
public class JobTest extends TestCase {

  private static class MockJobId extends IntegralAbstract implements JobId<MockJobId> {
    private final Integer mI;

    public MockJobId(final int i) {
      mI = i;
    }

    @Override
    public boolean validArguments(Result[] results) {
      return false;
    }

    @Override
    public boolean validResult(Result result) {
      return false;
    }

    @Override
    public boolean integrity() {
      return true;
    }

    @Override
    public int compareTo(MockJobId arg0) {
      return mI.compareTo(arg0.mI);
    }

    @Override
    public String toString() {
      return mI.toString();
    }

    @Override
    public boolean equals(Object arg0) {
      if (arg0 == null || !(arg0 instanceof MockJobId)) {
        return false;
      }
      return 0 == compareTo((MockJobId) arg0);
    }

    @Override
    public int hashCode() {
      return super.hashCode();
    }

    @Override
    public int time() {
      return mI;
    }

  }

  public void test() {
    final Job<MockJobId> job = new Job<MockJobId>(new MockJobId(42)) {
      @Override
      public Result run() {
        return null;
      }
    };
    assertEquals(42, job.id().mI.intValue());
    assertEquals("42", job.toString());
  }
}
