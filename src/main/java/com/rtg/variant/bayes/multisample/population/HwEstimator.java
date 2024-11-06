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

import java.util.List;

import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.HypothesisScores;
import com.rtg.variant.bayes.multisample.MultisampleJointScorer;
import com.rtg.variant.bayes.multisample.PriorContainer;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Records results of one iteration of EM algorithm for population estimators.
 * Assumes Hardy-Weinberg equilibrium and so needs only haploid counts.
 * Thus estimates fewer parameters than <code>DiploidEstimator</code> and should converge faster.
 */
public class HwEstimator implements Estimator {

  private static final double LAPLACE_CORRECTION = 1.0; //Double.parseDouble(System.getProperty("rtg.hwestimator.laplace-corr", "1.0"));


  final MultisampleJointScorer mScorer;

  /**
   * Create an estimator with a null joint scorer
   */
  public HwEstimator() {
    mScorer = new NullJointScorer();
  }

  HwEstimator(MultisampleJointScorer scorer) {
    mScorer = scorer;
  }

  @Override
  public <D extends Description, T extends HypothesesPrior<D>> EmResult<HypothesesPrior<D>> estimate(List<ModelInterface<?>> models, PriorContainer<T> priorContainer) {
    final HaploidDiploidHypotheses<T> hypotheses = priorContainer.getHypotheses();

    final HypothesisScores calls = mScorer.getBestScores(models, priorContainer);
    final int[] haploidCounts = new int[hypotheses.haploid().size()];

    int haploidTotal = 0;

    // initialise haploid counts to contain population allele counts when using population-priors. It is hard-to-impossible for complex calls.
    final DescriptionCounts dc = hypotheses.getDescriptionCounts();
    if (dc != null) {
      for (int i = 0; i < haploidCounts.length; ++i) {
        haploidCounts[i] += dc.getCount(i);
      }
      haploidTotal += dc.getTotalCount();
    }

    // TODO we should probably really only include calls from founders (and possibly half-founders), e.g. in order to exclude inherited de-novos from higher in the pedigree
    final Code code = hypotheses.diploid().code();
    for (int i = 0; i < models.size(); ++i) {
      final ModelInterface<?> model = models.get(i);
      final HypothesisScore call = calls.getScores()[i];
      if (call == null) { // E.g. for Female on Y chromosome or call from disagreeing hypotheses
        continue;
      }
      if (call.isDeNovo() == VariantSample.DeNovoStatus.IS_DE_NOVO) {
        continue;
      }
      final int best = call.hypothesis();
      if (model.haploid()) {
        haploidCounts[best]++;
        ++haploidTotal;
      } else {
        final int a = code.a(best);
        final int b = code.bc(best);
        haploidCounts[a]++;
        haploidCounts[b]++;
        haploidTotal += 2;
      }
    }
    final HaploidDiploidHypotheses<HypothesesPrior<D>> newHypothesis = computeNewPriors(hypotheses, haploidCounts, haploidTotal);
    return new EmResult<>(newHypothesis, calls, calls.getBs());
  }

  /**
   * Adjusts hypothesis priors according to Hardy-Weinberg principle from allele counts
   * @param hypotheses hypotheses containing input priors
   * @param haploidCounts counts for individual alleles, already including population level alleles
   * @param haploidTotal total count of alleles
   * @param <D> type of description
   * @param <T> type of hypotheses
   * @return list containing new hypothesis, first haploid, second diploid
   */
  public <D extends Description, T extends HypothesesPrior<D>> HaploidDiploidHypotheses<HypothesesPrior<D>> computeNewPriors(HaploidDiploidHypotheses<T> hypotheses, int[] haploidCounts, int haploidTotal) {
    final double laplaceCorrection = LAPLACE_CORRECTION;
    final PossibilityArithmetic arith = hypotheses.haploid().arithmetic();
    final T haploid = hypotheses.haploid();
    final T diploid = hypotheses.diploid();
    final Code code = diploid.code();

    //Compute new haploid priors.
    final double[] haploidProb = new double[haploidCounts.length];
    for (int i = 0; i < haploidCounts.length; ++i) {
      haploidProb[i] = (haploidCounts[i] + laplaceCorrection * arith.poss2Prob(haploid.p(i))) / (haploidTotal + laplaceCorrection);
    }
    final double[] haploidPoss = VariantUtils.prob2Poss(haploidProb, arith);
    final HypothesesPrior<D> newHaploid = new HypothesesPrior<>(hypotheses.haploid().description(), arith, haploidPoss, true, hypotheses.haploid().reference());

    //Compute new diploid priors.
    final int n = haploidTotal + 1;
    final double n2 = n * (double) n;
    final double[] diploidProb = new double[diploid.size()];
    for (int i = 0; i < diploidProb.length; ++i) {
      final int a = code.a(i);
      final int b = code.bc(i);

      final double haploidCorr = (a == b ? 1.0 : 2.0)
          * (haploidCounts[a] * (double) haploidCounts[b]
             + haploidCounts[a] * arith.poss2Prob(haploid.p(b))
             + haploidCounts[b] * arith.poss2Prob(haploid.p(a)));
      diploidProb[i] = (haploidCorr + arith.poss2Prob(diploid.p(i))) / n2;
    }
    final double[] diploidPoss = VariantUtils.prob2Poss(diploidProb, arith);
    final HypothesesPrior<D> newDiploid = new HypothesesPrior<>(hypotheses.diploid().description(), arith, diploidPoss, false, hypotheses.diploid().reference());

    return new HaploidDiploidHypotheses<>(hypotheses.none(), newHaploid, newDiploid, false, hypotheses.getDescriptionCounts());
  }

  @Override
  public String toString() {
    return "HwEstimator";
  }

}
