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
package com.rtg.variant.bayes.complex;

import com.rtg.util.integrity.Exam;

/**
 * Utilities for manipulating probability distributions.
 */
public final class IndelUtils {

  private IndelUtils() { } //prevent instantiation

  /**
   * Take the prob distribution whose first entry is 0.0 and create another normalized distribution
   * where the first entry is not and all the other entries have been reduced by a factor of scale.
   * @param prob probability distribution which is to be flattened.
   * @param scale by which to reduce the tail.
   * @return the flattened probability distribution.
   */
  static double[] flatten(final double[] prob, final double scale) {
    assert Exam.assertDistribution(prob);
    assert prob[0] == 0.0;
    assert 0.0 <= scale && scale <= 1.0 && !Double.isNaN(scale);
    final double[] res = new double[prob.length];
    res[0] = 1.0 - scale;
    double norm = 0.0;
    for (int i = 1; i < prob.length; ++i) {
      norm += prob[i];
    }
    for (int i = 1; i < prob.length; ++i) {
      res[i] = prob[i] * scale / norm;
    }
    assert Exam.assertDistribution(res);
    return res;
  }

  /**
   * Take the prob distribution and create another normalized distribution
   * where the first entry is 0.0 and all the other entries have been correspondingly scaled up.
   * @param prob probability distribution which is to be unflattened.
   * @return the unflattened probability distribution.
   */
  static double[] unFlatten(final double[] prob) {
    assert Exam.assertDistribution(prob);
    final double[] res = new double[prob.length];
    res[0] = 0.0;
    final double norm = 1.0 - prob[0];
    for (int i = 1; i < prob.length; ++i) {
      res[i] = prob[i] / norm;
    }
    assert Exam.assertDistribution(res);
    return res;
  }

}
