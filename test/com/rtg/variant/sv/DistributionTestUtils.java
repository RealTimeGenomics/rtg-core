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

package com.rtg.variant.sv;

import junit.framework.Assert;

/**
 */
public final class DistributionTestUtils {

  private DistributionTestUtils() { }

  static int peak(Distribution di) {
    double max = Double.NEGATIVE_INFINITY;
    int p = Integer.MIN_VALUE;
    for (int i = di.lo(); i < di.hi(); i++) {
      if (di.get(i) > max) {
        max = di.get(i);
        p = i;
      }
    }
    return p;
  }

  static void checkMin(final SamCounts sa, final Distribution distr, final int expMini, final double expMin) {
    final Signal sig = distr.getSignalLn(sa, "");
    double min = Double.POSITIVE_INFINITY;
    int mini = -1;
    for (int i = 0; i < 20; i++) {
      final double v = sig.value(i);
      if (v < min) {
        min = v;
        mini = i;
      }
      //System.err.println(i + " " + v);
    }
    Assert.assertEquals(expMini, mini);
    Assert.assertEquals(expMin, min, 0.00001);
  }


}
