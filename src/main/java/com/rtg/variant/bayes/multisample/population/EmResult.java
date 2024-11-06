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

package com.rtg.variant.bayes.multisample.population;

import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.HypothesisScores;
import com.rtg.variant.bayes.multisample.PriorContainer;
import com.rtg.variant.bayes.multisample.forwardbackward.BContainer;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * Result of one iteration of EM algorithm.
 * @param <T> type of hypotheses
 */
public class EmResult<T extends HypothesesPrior<?>> {

  private final HypothesisScores mCalls;

  private final PriorContainer<T> mPriorContainer;


  /**
   * @param hypotheses the haploid diploid hypotheses
   * @param calls made by the models.
   * @param bs B values for all samples
   */
  public EmResult(HaploidDiploidHypotheses<T> hypotheses, HypothesisScores calls, BContainer[] bs) {
    mPriorContainer = new PriorContainer<>(hypotheses, bs);
    mCalls = calls;
  }

  PriorContainer<T> getPriorContainer() {
    return mPriorContainer;
  }

  int difference(final HypothesisScores calls) {
    return difference(calls.getScores(), mCalls.getScores());
  }

  static int difference(final HypothesisScore[] calls0, final HypothesisScore[] calls1) {
    assert calls0.length == calls1.length;
    int count = 0;
    for (int i = 0; i < calls0.length; ++i) {
      if (calls0[i] == null) { // E.g. female on Y chromosome
        if (calls1[i] != null) { // Can't suddenly change across em iterations.
          throw new RuntimeException();
        }
      } else if (calls0[i].hypothesis() != calls1[i].hypothesis()) {
        ++count;
      }
    }
    return count;
  }

  /**
   * Get calls.
   * @return Returns the calls.
   */
  HypothesisScores calls() {
    return mCalls;
  }

}
