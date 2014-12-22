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

package com.rtg.variant.bayes;

import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Provides a probability distribution over a (haploid) set of hypotheses
 * that provide evidence for Bayesian estimation.
 * They will sum to approximately 1.0 and are in standard double format.
 * (At this level we are not using <code>PossibilityArithmetic</code>).
 */
public abstract class Evidence extends IntegralAbstract implements EvidenceInterface {
  private final Description mDescription;

  private final double mMapError;

  /**
   * @param description of the haploid hypotheses.
   * @param mapError probability that the read doesn't map to this position.
   */
  protected Evidence(final Description description, final double mapError) {
    mDescription = description;
    mMapError = mapError;
  }

  /**
   * @return description of the haploid hypotheses.
   */
  protected Description description() {
    return mDescription;
  }

  @Override
  public double mapError() {
    return mMapError;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    double sum = 0.0;
    for (int i = 0; i < mDescription.size(); i++) {
      final double p = probability(i);
      Exam.assertTrue(0.0 <= p && p <= 1.0 && !Double.isNaN(p));
      sum += p;
    }
    Exam.assertEquals(1.0, sum, 0.000001);
    return true;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append("Evidence ");
    sb.append(isForward() ? ">" : "<");
    sb.append(" pe=").append(Utils.realFormat(pe(), 3));
    sb.append(" error=").append(Utils.realFormat(error(), 3));
    sb.append(" notMap=").append(Utils.realFormat(mapError(), 3));
    sb.append(LS);
    final int pad = mDescription.maxLength();
    for (int i = 0; i < mDescription.size(); i++) {
      sb.append(" ");
      sb.append(StringUtils.padLeft(mDescription.name(i), pad));
      sb.append(" ");
      final double probability = Math.log(probability(i));
      sb.append(StringUtils.padLeft(Utils.realFormat(probability, 3), 6));
      sb.append(LS);
    }
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mDescription);
    return true;
  }

}
