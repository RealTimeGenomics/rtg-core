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

  private JobType(int numberArgs) {
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
    for (int i = 0; i < mNumberArgs; i++) {
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
