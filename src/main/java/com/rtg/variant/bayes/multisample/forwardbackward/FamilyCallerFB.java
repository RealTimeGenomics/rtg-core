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
package com.rtg.variant.bayes.multisample.forwardbackward;


import java.util.Arrays;
import java.util.List;

import com.rtg.relation.Family;
import com.rtg.util.MathUtils;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Factor;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCaller;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.HypothesisScores;
import com.rtg.variant.bayes.multisample.MultisampleJointScorer;
import com.rtg.variant.bayes.multisample.PriorContainer;
import com.rtg.variant.bayes.multisample.family.MendelianAlleleProbabilityFactory;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * Generates the family comparison. Assumes the output will come in with
 * father and mother bayesians in first two array positions.
 */
public class FamilyCallerFB extends AbstractMultisampleCaller implements MultisampleJointScorer {

  final VariantParams mParams;
  private final Family[] mFamilies;

  /**
   * @param params genome priors
   * @param families the family structures
   */
  public FamilyCallerFB(VariantParams params, Family... families) {
    mParams = params;
    mFamilies = families;
  }

  // Put a new score into the current scores array, dealing with cases where we have already called the individual
  static void setScore(HypothesisScore[] scores, int index, HypothesisScore newScore) {
    if (scores[index] == null) {
      scores[index] = newScore;
    } else {
      // FB in theory can't produce disagreeing hypotheses. We'll now be confirming this.
      if (scores[index].hypothesis() != newScore.hypothesis()) {
        throw new UnsupportedOperationException("disagreeing calls in FB.");
      }
    }
  }

  /**
   * Generate best scores for a family
   *
   * @param <D> the type of the description.
   * @param <T> the type of the hypotheses prior.
   * @param models input models to call from
   * @param priorContainer hypotheses and <code>B</code>
   * @return the scores.
   */
  @Override
  public <D extends Description, T extends HypothesesPrior<D>> HypothesisScores getBestScores(List<ModelInterface<?>> models, PriorContainer<T> priorContainer) {
    boolean isInteresting = false;
    final HypothesisScore[] scores = new HypothesisScore[models.size()];
    final int numberIterations = 2;
    double sumNips = 0;
    for (int oneloop = 0; oneloop < numberIterations; ++oneloop) {
      final FamilyPosteriorFB[] fps = new FamilyPosteriorFB[mFamilies.length];
      final Factor<?>[] as = CommonFormulas.initialA(models, priorContainer.getHypotheses());
      final Factor<?>[] asMendel = CommonFormulas.initialA(models, priorContainer.getHypotheses());
      final Factor<?>[] asDenovo = CommonFormulas.initialA(models, priorContainer.getHypotheses());
      final FamilyPosteriorFB[] fpsMendel = new FamilyPosteriorFB[mFamilies.length];
      final FamilyPosteriorFB[] fpsDenovo = new FamilyPosteriorFB[mFamilies.length];
      final BContainer[] bs = priorContainer.getBs();
      assert bs != null;
      assert models.size() == bs.length;
      final MendelianAlleleProbabilityFactory fact = MendelianAlleleProbabilityFactory.COMBINED;
      for (int familyIndex = mFamilies.length - 1; familyIndex >= 0; --familyIndex) {
        final Family family = mFamilies[familyIndex];
        final FamilyPosteriorFB fp = new FamilyPosteriorFB(family, mParams.genomePriors(), models, priorContainer.getHypotheses(), fact);
        fps[familyIndex] = fp;
        final int[] ids = family.getSampleIds();
        final Factor<?>[] newBs = fp.computeParentBs(as, bs);
        bs[ids[Family.FATHER_INDEX]].setB(family.getFatherFamilyId(), newBs[Family.FATHER_INDEX]);
        bs[ids[Family.MOTHER_INDEX]].setB(family.getMotherFamilyId(), newBs[Family.MOTHER_INDEX]);
      }

      for (int familyIndex = 0; familyIndex < mFamilies.length; ++familyIndex) {
        final Family family = mFamilies[familyIndex];
        final FamilyPosteriorFB fp = fps[familyIndex]; // new FamilyPosteriorFB(family, mParams.genomePriors(), models, priorContainer.getHypotheses(), mParams.denovoPrior());
        final int[] ids = family.getSampleIds();

        final Factor<?>[] newAs = fp.computeChildAs(as, bs);

        //denovo handling
        final FamilyPosteriorFB fpMendel = new FamilyPosteriorFB(family, mParams.genomePriors(), models, priorContainer.getHypotheses(), MendelianAlleleProbabilityFactory.MENDELIAN);
        final Factor<?>[] asMendNew = fpMendel.computeChildAs(as, bs);
        final FamilyPosteriorFB fpDenovo = new FamilyPosteriorFB(family, mParams.genomePriors(), models, priorContainer.getHypotheses(), MendelianAlleleProbabilityFactory.DENOVO);
        final Factor<?>[] asDenovoNew = fpDenovo.computeChildAs(as, bs);
        fpsMendel[familyIndex] = fpMendel;
        fpsDenovo[familyIndex] = fpDenovo;
        //

        // todo, this implementation does not have any correction for contrary evidence

        for (int c = 0; c < family.numChildren(); ++c) {
          as[ids[Family.FIRST_CHILD_INDEX + c]] = newAs[c]; //TODO should we update in the method itself?
        }

        //more denovo handling
        for (int c = 0; c < family.numChildren(); ++c) {
          asMendel[ids[Family.FIRST_CHILD_INDEX + c]] = asMendNew[c];
        }
        for (int c = 0; c < family.numChildren(); ++c) {
          asDenovo[ids[Family.FIRST_CHILD_INDEX + c]] = asDenovoNew[c];
        }
      }

      if (oneloop + 1 >= numberIterations) {
//        System.err.print("FCFB :");
        for (int familyIndex = 0; familyIndex < mFamilies.length; ++familyIndex) {
          final Family family = mFamilies[familyIndex];
          final FamilyPosteriorFB fp = fps[familyIndex];
          final int[] ids = family.getSampleIds();
          //even more denovo stuff
          final FamilyPosteriorFB mendelFamily = fpsMendel[familyIndex];
          final FamilyPosteriorFB denovoFamily = fpsDenovo[familyIndex];
          mendelFamily.computeMarginals(asMendel, bs);
          denovoFamily.computeMarginals(asDenovo, bs);
          final double[] childDenovoPosteriors = new double[family.numChildren()];
          for (int i = 0; i < family.numChildren(); ++i) {
            childDenovoPosteriors[i] = denovoFamily.childMarginalSum(i) - mendelFamily.childMarginalSum(i);
          }

          fp.addMarginals(mendelFamily);
          fp.addMarginals(denovoFamily);
          fp.findBest(childDenovoPosteriors);
          //
          setScores(scores, family, fp, ids);
          if (fp.isInteresting()) {
            isInteresting = true;
          }
          final double nip = fp.getNonIdentityPosterior();
          sumNips += MathUtils.logExpPlus1(nip); // Math.log(Math.exp(nip) + 1), nip = -ip for ref= calls
        }
//        System.err.println(" -> " + sumNips);
      }
    }
    final double qual = MathUtils.logExpMinus1(sumNips); // Math.log(Math.exp(sumNips) - 1.0);
    return new HypothesisScores(scores, isInteresting, qual, priorContainer.getBs());
  }

  private void setScores(HypothesisScore[] scores, Family family, FamilyPosteriorFB fp, int[] ids) {
    setScore(scores, ids[Family.FATHER_INDEX], fp.bestFather());
    setScore(scores, ids[Family.MOTHER_INDEX], fp.bestMother());
    for (int i = 0; i < family.numChildren(); ++i) {
      setScore(scores, ids[Family.FIRST_CHILD_INDEX + i], fp.bestChild(i));
    }
  }

  @Override
  public BContainer[] makeInitialBs(List<ModelInterface<?>> models) {
    final int[] bsizes = new int[models.size()];
    Arrays.fill(bsizes, 1);
    for (Family f : mFamilies) {
//      if (f.getFatherDistinctMates() > 1 || f.getMotherDistinctMates() > 1) {
//        throw new UnsupportedOperationException("detected bsize > 1");
//      }
      bsizes[f.getSampleIds()[Family.FATHER_INDEX]] = f.getFatherDistinctMates();
      bsizes[f.getSampleIds()[Family.MOTHER_INDEX]] = f.getMotherDistinctMates();
    }
    return CommonFormulas.initialB(models, bsizes);
  }

  @Override
  protected <D extends Description, T extends HypothesesPrior<D>> ComparisonResult makeSamples(List<ModelInterface<?>> models, HaploidDiploidHypotheses<T> hypotheses) {
    //this isn't used in the population caller context, so should never be called
    throw new UnsupportedOperationException();
  }

  @Override
  protected VariantParams getParams() {
    return mParams;
  }
}
