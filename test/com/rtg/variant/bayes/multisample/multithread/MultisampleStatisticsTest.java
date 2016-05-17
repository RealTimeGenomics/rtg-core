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

import junit.framework.TestCase;

/**
 */
public class MultisampleStatisticsTest extends TestCase {

  public void test() {
    final MultisampleStatistics st = new MultisampleStatistics();
    st.increment(new JobIdMultisample(5, 0, JobType.INCR), 42);
    assertEquals("[Timer INCR      0.00  count 1       0.00 bytes read 0, Timer DANGLING empty, Timer COMPLEX empty, Timer FLUSH empty, Timer FILTER empty, Timer BED empty, Timer OUT empty]"
        , st.toString());
  }
}
