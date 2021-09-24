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

import java.util.Arrays;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Objects returned by <code>Job</code>s.
 */
public class Result extends IntegralAbstract {

  private final Object[] mResults;

  /**
   * @param results from a job.
   */
  public Result(final Object... results) {
    mResults = results;
  }

  /**
   * @param index selects the result.
   * @return an individual result.
   */
  public Object result(final int index) {
    return mResults[index];
  }

  /**
   * @return the number of individual results.
   */
  public int length() {
    return mResults.length;
  }

  @Override
  public String toString() {
    return Arrays.toString(mResults);
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mResults);
    return true;
  }


}
