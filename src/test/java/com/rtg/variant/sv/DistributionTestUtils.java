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

package com.rtg.variant.sv;

import org.junit.Assert;

/**
 */
public final class DistributionTestUtils {

  private DistributionTestUtils() { }

  static int peak(Distribution di) {
    double max = Double.NEGATIVE_INFINITY;
    int p = Integer.MIN_VALUE;
    for (int i = di.lo(); i < di.hi(); ++i) {
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
    for (int i = 0; i < 20; ++i) {
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
