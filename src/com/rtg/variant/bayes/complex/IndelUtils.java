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
