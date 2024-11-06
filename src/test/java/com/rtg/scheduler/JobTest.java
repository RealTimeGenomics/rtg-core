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

package com.rtg.scheduler;

import com.rtg.util.integrity.IntegralAbstract;

import junit.framework.TestCase;

/**
 */
public class JobTest extends TestCase {

  private static class MockJobId extends IntegralAbstract implements JobId<MockJobId> {
    private final Integer mI;

    MockJobId(final int i) {
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
