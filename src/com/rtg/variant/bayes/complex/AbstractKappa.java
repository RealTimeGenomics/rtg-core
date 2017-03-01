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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
@TestClass(value = {"com.rtg.variant.bayes.complex.KappaImplementationTest", "com.rtg.variant.bayes.complex.KappaMakeTester"})
public abstract class AbstractKappa extends IntegralAbstract {

  /**
   * The proper way of getting an {@link AbstractKappa} constructed.
   * @param insEventRate event rate for insertions.
   * @param insDistribution distribution of insertion lengths. (0 length assumed 0.0)
   * @param delEventRate event rate for deletions.
   * @param delDistribution distribution of deletion lengths. (0 length assumed 0.0).
   * @param indelDecay rate at which indel priors decrease beyond the end of the explicit distributions.
   * @return an implementation of <code>AbstractKappa</code>
   */
  public static AbstractKappa makeKappa(final double insEventRate, final double[] insDistribution, final double delEventRate, final double[] delDistribution, final double indelDecay) {
    return new KappaMemo(new KappaImplementation(insEventRate, insDistribution, delEventRate, delDistribution, indelDecay));
  }

  /**
   * Return the two parameter version of kappa (see theory).
   * @param m length of evidence.
   * @param l length of hypothesis.
   * @return the value of kappa.
   */
  public double kappa(final int m, final int l) {
    if (m < 0) {
      throw new IndexOutOfBoundsException(String.valueOf(m));
    }
    if (l < 0) {
      throw new IndexOutOfBoundsException(String.valueOf(l));
    }
    final double pi = pi(m - l);
    final double corr = piSum(-l);
    //System.err.println("m=" + m + " l=" + l + " pi=" + pi + " piSum=" + corr);
    final double kappa = pi / corr;
    assert kappa >= 0.0 && kappa <= 1.0 && !Double.isNaN(kappa);
    return kappa;
  }


  /**
   * Return pi (see theory).
   * @param i position in pi (often computed as m - l)
   * @return the value of pi.
   */
  abstract double pi(int i);

  /**
   * Return the cumulative sum of pi (from l to positive infinity).
   * @param i position in <code>piSum</code>.
   * @return the value of sum from i (inclusive) to infinity..
   */
  abstract double piSum(int i);

}
