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


import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import com.rtg.util.ChiSquared;
import com.rtg.util.Utils;

/**
 * Used to generate distributions.
 */
public final class DistributionUtils {

  private DistributionUtils() { }

  /**
   * A cumulative distribution over a normal distribution.
   * @param lo low index in distribution.
   * @param hi high index in distribution.
   * @param rate the maximum rate multiplied by the distribution.
   * @param mean of normal distribution.
   * @param stdDev standard deviation of normal distribution.
   * @param reverse if true return the reverse distribution.
   * @return the distribution.
   */
  static Distribution distribution(final int lo, final int hi, final double rate, final double mean, final double stdDev, boolean reverse) {
    if (reverse) {
      return distribution2(lo, hi, rate, mean + 1, stdDev);
    } else {
      return distribution1(lo, hi, rate, mean, stdDev);
    }
  }

  /**
   * A cumulative distribution over a normal distribution.
   * constant to the left and decreasing to 0.0 on the right.
   * @param lo low index in distribution.
   * @param hi high index in distribution.
   * @param rate the maximum rate multiplied by the distribution.
   * @param mean of normal distribution.
   * @param stdDev standard deviation of normal distribution.
   * @return the distribution.
   */
  static Distribution distribution1(final int lo, final int hi, final double rate, final double mean, final double stdDev) {
    final int diameter = hi - lo;
    final double[] distr = new double[diameter];
    for (int i = 0; i < diameter; ++i) {
      final int j = i + lo;
      final double diff = mean - j;
      final double x = diff / stdDev;
      distr[i] = rate * ChiSquared.normal(x);
    }
    return new DistributionArray(lo, distr);
  }

  /**
   * A cumulative distribution over a normal distribution.
   * constant to the right and decreasing to 0.0 on the left.
   * @param lo low index in distribution.
   * @param hi high index in distribution.
   * @param rate the maximum rate multiplied by the distribution.
   * @param mean of normal distribution.
   * @param stdDev standard deviation of normal distribution.
   * @return the distribution.
   */
  static Distribution distribution2(final int lo, final int hi, final double rate, final double mean, final double stdDev) {
    final int diameter = hi - lo;
    final double[] distr = new double[diameter];
    for (int i = 0; i < diameter; ++i) {
      final int j = i + lo;
      final double diff = mean - j;
      final double x = -diff / stdDev;
      distr[i] = rate * ChiSquared.normal(x);
    }
    return new DistributionArray(lo, distr);
  }

  /**
   * Subtract a distribution from a constant.
   * @param a constant being subtracted from (minuend).
   * @param b distribution being subtracted (subtrahend).
   * @return difference.
   */
  static Distribution subtract(double a, Distribution b) {
    final int lo = b.lo();
    final int diameter = b.hi() - lo;
    final double[] distr = new double[diameter];
    for (int i = 0; i < diameter; ++i) {
      final double d = a - b.get(i + lo);
      final double e;
      if (d < 0.0) {
        if (d > -0.000001) {
          //allow for rounding issues
          e = 0.0;
        } else {
          throw new RuntimeException("i=" + i + " a=" + Utils.realFormat(a, 7) + " b=" + Utils.realFormat(b.get(i + lo), 7) + " d=" + Utils.realFormat(d, 7));
        }
      } else {
        e = d;
      }
      distr[i] = e;
    }
    return new DistributionArray(lo, distr);
  }

  /**
   * Add a constant to a distribution.
   * @param a summand distribution.
   * @param b constant being added.
   * @return sum.
   */
  static Distribution add(Distribution a, double b) {
    final int lo = a.lo();
    final int hi = a.hi();
    final int diameter = hi - lo;
    final double[] distr = new double[diameter];
    for (int i = lo; i < hi; ++i) {
      distr[i - lo] = a.get(i) + b;
    }
    return new DistributionArray(lo, distr);
  }

  /**
   * Pointwise sum of two distributions.
   * @param a first distribution.
   * @param b second distribution.
   * @return sum distribution.
   */
  static Distribution add(Distribution a, Distribution b) {
    final int lo = a.lo();
    final int hi = a.hi();
    final int diameter = hi - lo;
    assert lo == b.lo();
    assert hi == b.hi() : hi + ":" + b.hi();
    final double[] distr = new double[diameter];
    for (int i = lo; i < hi; ++i) {
      distr[i - lo] = a.get(i) + b.get(i);
    }
    return new DistributionArray(lo, distr);
  }

  /**
   * Pointwise product of two distributions.
   * @param a first distribution.
   * @param b second distribution.
   * @return product distribution.
   */
  static Distribution multiply(Distribution a, Distribution b) {
    final int lo = a.lo();
    final int hi = a.hi();
    final int diameter = hi - lo;
    assert lo == b.lo();
    assert hi == b.hi();
    final double[] distr = new double[diameter];
    for (int i = lo; i < hi; ++i) {
      distr[i - lo] = a.get(i) * b.get(i);
    }
    return new DistributionArray(lo, distr);
  }

  static void plot(File outFile, String id, Distribution distr) {
    try {
      try (PrintStream ps = new PrintStream(new FileOutputStream(outFile))) {
        plot(id, distr, ps);
      }
    } catch (final IOException e) {
      throw new RuntimeException("Could not generate plot for " + id, e);
    }
  }

  static void plot(final String id, final Distribution distr, final PrintStream out) {
    out.println();
    out.println("#" + id);
    out.print(distr.dump());
  }
}
