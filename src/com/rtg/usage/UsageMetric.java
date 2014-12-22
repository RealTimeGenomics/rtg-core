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
package com.rtg.usage;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class UsageMetric extends IntegralAbstract {
  private long mMetric = 0;

  /**
   * @return the value of the metric counter.
   */
  public synchronized long getMetric() {
    return mMetric;
  }

  /**
   * This may only be called once.
   * @param metric the value to set the metric to.
   */
  public synchronized void setMetric(long metric) {
    assert metric >= 0;
    mMetric = metric;
  }

  /**
   * @param metric value to increment the metric by.
   */
  public synchronized void incrementMetric(long metric) {
    assert metric >= 0;
    mMetric += metric;
  }

  @Override
  public  synchronized boolean integrity() {
    Exam.assertTrue(mMetric >= -1);
    return true;
  }
}
