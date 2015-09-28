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
    for (int oneloop = 0; oneloop < numberIterations; oneloop++) {
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
      for (int familyIndex = mFamilies.length - 1; familyIndex >= 0; familyIndex--) {
        final Family family = mFamilies[familyIndex];
        final FamilyPosteriorFB fp = new FamilyPosteriorFB(family, mParams.genomePriors(), models, priorContainer.getHypotheses(), fact);
        fps[familyIndex] = fp;
        final int[] ids = family.getSampleIds();
        final Factor<?>[] newBs = fp.computeParentBs(as, bs);
        bs[ids[Family.FATHER_INDEX]].setB(family.getFatherFamilyId(), newBs[Family.FATHER_INDEX]);
        bs[ids[Family.MOTHER_INDEX]].setB(family.getMotherFamilyId(), newBs[Family.MOTHER_INDEX]);
      }

      for (int familyIndex = 0; familyIndex < mFamilies.length; familyIndex++) {
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

        for (int c = 0; c < family.numChildren(); c++) {
          as[ids[Family.FIRST_CHILD_INDEX + c]] = newAs[c]; //TODO should we update in the method itself?
        }

        //more denovo handling
        for (int c = 0; c < family.numChildren(); c++) {
          asMendel[ids[Family.FIRST_CHILD_INDEX + c]] = asMendNew[c];
        }
        for (int c = 0; c < family.numChildren(); c++) {
          asDenovo[ids[Family.FIRST_CHILD_INDEX + c]] = asDenovoNew[c];
        }
      }

      if (oneloop + 1 >= numberIterations) {
//        System.err.print("FCFB :");
        for (int familyIndex = 0; familyIndex < mFamilies.length; familyIndex++) {
          final Family family = mFamilies[familyIndex];
          final FamilyPosteriorFB fp = fps[familyIndex];
          final int[] ids = family.getSampleIds();
          //even more denovo stuff
          final FamilyPosteriorFB mendelFamily = fpsMendel[familyIndex];
          final FamilyPosteriorFB denovoFamily = fpsDenovo[familyIndex];
          mendelFamily.computeMarginals(asMendel, bs);
          denovoFamily.computeMarginals(asDenovo, bs);
          final double[] childDenovoPosteriors = new double[family.numChildren()];
          for (int i = 0; i < family.numChildren(); i++) {
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
    for (int i = 0; i < family.numChildren(); i++) {
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
