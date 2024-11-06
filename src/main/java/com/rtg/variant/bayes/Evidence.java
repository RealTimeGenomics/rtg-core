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
    for (int i = 0; i < mDescription.size(); ++i) {
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
    for (int i = 0; i < mDescription.size(); ++i) {
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
