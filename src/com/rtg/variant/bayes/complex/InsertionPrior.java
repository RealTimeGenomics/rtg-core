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

import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.GenomePriorParams;

/**
 */
public class InsertionPrior extends IntegralAbstract {

  private static AbstractKappa flatKappa(final double[] distr, double indelDecay) {
    assert Exam.assertDistribution(distr);
    final double zRate = (1.0 - distr[0]) / (2 - distr[0]);
    final double[] zetaRemainder = IndelUtils.unFlatten(distr);
    return AbstractKappa.makeKappa(
        zRate,
        zetaRemainder,
        zRate,
        zetaRemainder,
        indelDecay
    );
  }

  private final double[] mZetaPrime;

  private final double[] mZeta;

  private final AbstractKappa mHeteroPriorKappa;

  private final AbstractKappa mHomoPriorKappa;

  /**
   * @param params parameters
   */
  public InsertionPrior(final GenomePriorParams params) {
    //zeta is a flattened distribution (mZeta[0] > 0) for use with heterozygous priors
    mZeta = params.genomeIndelDistribution();
    assert Exam.assertDistribution(mZeta);
    mHeteroPriorKappa = flatKappa(mZeta, params.genomeIndelLengthDecay());

    final double genomeIndelEventRate = params.genomeIndelEventRate();
    final double[] zetaUnFlat = IndelUtils.unFlatten(mZeta);
    //System.err.println("zeta=" + Arrays.toString(mZeta));
    //System.err.println("zetaUnFlat=" + Arrays.toString(zetaUnFlat));
    mZetaPrime = IndelUtils.flatten(zetaUnFlat, genomeIndelEventRate);
    //System.err.println("zetaPrime=" + Arrays.toString(mZetaPrime));
    assert Exam.assertDistribution(mZetaPrime);
    mHomoPriorKappa = flatKappa(mZetaPrime, params.genomeIndelLengthDecay());

    assert integrity();
  }

  /**
   * Get the kappa to be used for heterozygous prior calculations.
   * @return the kappa to be used for heterozygous prior calculations.
   */
  public AbstractKappa heteroPriorKappa() {
    return mHeteroPriorKappa;
  }

  /**
   * Get the kappa to be used for homozygous prior calculations.
   * @return the kappa to be used for homozygous prior calculations.
   */
  public AbstractKappa homoPriorKappa() {
    return mHomoPriorKappa;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("Insertion prior").append(LS);
    sb.append("zeta          ").append(Utils.realFormat(mZeta, 6)).append(LS);
    sb.append("zeta`         ").append(Utils.realFormat(mZetaPrime, 6)).append(LS);
    sb.append(LS);
  }

  @Override
  public boolean integrity() {
    //System.err.println(toString());
    //zeta`
    Exam.assertDistribution(mZetaPrime);
    //zeta
    Exam.assertDistribution(mZeta);
    return true;
  }
}
