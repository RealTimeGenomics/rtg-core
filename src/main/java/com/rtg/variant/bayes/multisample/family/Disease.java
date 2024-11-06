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
package com.rtg.variant.bayes.multisample.family;

import com.rtg.relation.Family;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Hypotheses;

/**
 * Disease consistency checking. See <code>Q</code> function in <code>multiScoring.pdf</code>
 */
public class Disease {
  /** A prediction that the disease is not explained by the current call. */
  private final Family mFamily;
  private final Hypotheses<?> mIndividualHypotheses;
  private final int[] mCurrentHypotheses;

  /**
   * @param family describes names of family and disease status.
   * @param individualHypotheses describes the hypotheses associated with each individual family member (NB: NOT a disease hypothesis).
   * @param currentHypotheses the current hypotheses associated with each family member as described by <code>individualHypotheses</code>.
   */
  //TODO add recessive vs dominant flag
  Disease(final Family family, final Hypotheses<?> individualHypotheses, final int... currentHypotheses) {
    if (!family.isOneParentDiseased()) {
      throw new UnsupportedOperationException();
    }
    mFamily = family;
    mIndividualHypotheses = individualHypotheses;
    mCurrentHypotheses = currentHypotheses;
  }

  /**
   * Find an explanation. Because of the assumptions about families (one parent diseased, another not and 1 or more children)
   * there will either be one explanation or none.
   * @param ref reference in code space of individual hypotheses.
   * @return a disease explanation in the code space of the diseases (see <code>HypothesesDisease</code>)
   */
  public int explanation(final int ref) {
    for (int i = 0; i < mIndividualHypotheses.description().size(); ++i) {
      if (i != ref && q(i + 1)) {
        return i + 1;
      }
    }
    return 0;
  }

  /**
   * @param diseaseHyp haploid disease hypothesis that starts with NONE and includes rest of nucleotides (see <code>HypothesesDisease</code>).
   * @return true iff <code>diseaseHyp</code> explains all the disease status.
   */
  private boolean q(final int diseaseHyp) {
    //iterate over the genomes
    for (int i = 0; i < mCurrentHypotheses.length; ++i) {
      if (!q(diseaseHyp, i)) {
        return false;
      }
    }
    return true;
  }

  /**
   * @param diseaseHyp haploid disease hypothesis that starts with NONE and includes rest of nucleotides (see <code>HypothesesDisease</code>).
   * @param genome selects family member with father starting at 0.
   * @return true iff the disease hypothesis is consistent with the disease status of this family member (aka genome). Assumes disease is dominant.
   */
  private boolean q(final int diseaseHyp, final int genome) {
    final int hyp = mCurrentHypotheses[genome];
    final boolean contains;
    final Code code = mIndividualHypotheses.code();
    final int dis = diseaseHyp - 1;
    //TODO implement recessive version controlled by flag.
    if (code.homozygous(hyp)) {
      contains = code.a(hyp) == dis;
    } else {
      contains = code.a(hyp) == dis || code.b(hyp) == dis;
    }
    return contains == mFamily.isDiseased(genome);
  }
}
