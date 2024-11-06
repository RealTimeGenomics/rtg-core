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

import com.rtg.scheduler.Result;

/**
 * Execution phases within one chunk.
 * These are dummies for the sake of example but I have included comments as if they were the real thing.
 */
public enum JobType {
  /** Increment model counts and do initial calling. */
  INCR(0),
  /** Tidy dangling regions at the ends of chunks. */
  DANGLING(2),
  /** Complex calling. */
  COMPLEX(1),
  /** Merge results from initial calling and complex calling. */
  MERGE(2),
  /** Write regions to cx.bed file. */
  BED(2),
  /** Write called SNPs. */
  OUT(2);

  private final int mNumberArgs;

  JobType(int numberArgs) {
    mNumberArgs = numberArgs;
  }

  /**
   * @param results to be validated.
   * @return true iff the number and types of the results are correct as arguments for a job with this identifier.
   */
  public boolean  validArguments(Result[] results) {
    if (results.length != mNumberArgs) {
      return false;
    }
    for (int i = 0; i < mNumberArgs; ++i) {
      if (!validResult(results[i])) {
        return false;
      }
    }
    return true;
  }

  /**
   * @param result to be validated.
   * @return true iff the number and types of the objects in result are correct as the result for a job with this identifier.
   */
  public boolean validResult(Result result) {
    //System.err.println("validResult " + this + ":" + result);
    return result == null || result.length() == 1 && result.result(0) instanceof String;
  }

}
