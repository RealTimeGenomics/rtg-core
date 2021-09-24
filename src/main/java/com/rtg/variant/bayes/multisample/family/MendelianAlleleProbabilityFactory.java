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
package com.rtg.variant.bayes.multisample.family;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.reference.Ploidy;
import com.rtg.util.Pair;
import com.rtg.variant.bayes.Code;

/**
 */
@TestClass(value = {"com.rtg.variant.bayes.multisample.family.FastFamilyPosteriorTest", "com.rtg.variant.bayes.multisample.forwardbackward.FamilyCallerFBTest"})
public abstract class MendelianAlleleProbabilityFactory {

  /** Probability tables for doing calls including mendelian and de novo probabilities */
  public static final MendelianAlleleProbabilityFactory COMBINED = new Combined();
  /** Probability tables for doing calls including just mendelian probabilities */
  public static final MendelianAlleleProbabilityFactory MENDELIAN = new Mendelian();
  /** Probability tables for doing calls including just de novo probabilities */
  public static final MendelianAlleleProbabilityFactory DENOVO = new Denovo();


  static final Pair<MendelianAlleleProbability, MendelianAlleleProbability> NONE_HAPLOID_NONE = new Pair<>(MendelianAlleleProbabilityNHN.SINGLETON_NH, MendelianAlleleProbabilityNHNDeNovo.SINGLETON_NH);
  static final Pair<MendelianAlleleProbability, MendelianAlleleProbability> NONE_HAPLOID_HAPLOID = new Pair<>(MendelianAlleleProbabilityNHH.SINGLETON_NH, MendelianAlleleProbabilityNHHDeNovo.SINGLETON_NH);
  static final Pair<MendelianAlleleProbability, MendelianAlleleProbability> HAPLOID_NONE_NONE = new Pair<>(MendelianAlleleProbabilityNHN.SINGLETON_HN, MendelianAlleleProbabilityNHNDeNovo.SINGLETON_HN);
  static final Pair<MendelianAlleleProbability, MendelianAlleleProbability> HAPLOID_NONE_HAPLOID = new Pair<>(MendelianAlleleProbabilityNHH.SINGLETON_HN, MendelianAlleleProbabilityNHHDeNovo.SINGLETON_HN);
  static final Pair<MendelianAlleleProbability, MendelianAlleleProbability> HAPLOID_HAPLOID_HAPLOID = new Pair<>(MendelianAlleleProbabilityHHH.SINGLETON_HHH, MendelianAlleleProbabilityHHHDeNovo.SINGLETON_HHH);
  static final Pair<MendelianAlleleProbability, MendelianAlleleProbability> HAPLOID_HAPLOID_DIPLOID = new Pair<>(MendelianAlleleProbabilityHHD.SINGLETON_HHD, MendelianAlleleProbabilityHHDDeNovo.SINGLETON_HHD);
  static final Pair<MendelianAlleleProbability, MendelianAlleleProbability> HAPLOID_DIPLOID_HAPLOID = new Pair<>(MendelianAlleleProbabilityHDH.SINGLETON_HD, MendelianAlleleProbabilityHDHDeNovo.SINGLETON_HD);
  static final Pair<MendelianAlleleProbability, MendelianAlleleProbability> HAPLOID_DIPLOID_DIPLOID = new Pair<>(MendelianAlleleProbabilityHDD.SINGLETON_HD, MendelianAlleleProbabilityHDDDeNovo.SINGLETON_HD);
  static final Pair<MendelianAlleleProbability, MendelianAlleleProbability> DIPLOID_HAPLOID_HAPLOID = new Pair<>(MendelianAlleleProbabilityHDH.SINGLETON_DH, MendelianAlleleProbabilityHDHDeNovo.SINGLETON_DH);
  static final Pair<MendelianAlleleProbability, MendelianAlleleProbability> DIPLOID_HAPLOID_DIPLOID = new Pair<>(MendelianAlleleProbabilityHDD.SINGLETON_DH, MendelianAlleleProbabilityHDDDeNovo.SINGLETON_DH);
  static final Pair<MendelianAlleleProbability, MendelianAlleleProbability> DIPLOID_DIPLOID_DIPLOID = new Pair<>(MendelianAlleleProbabilityDiploid.SINGLETON, MendelianAlleleProbabilityDiploidDeNovo.SINGLETON);


  private static Pair<MendelianAlleleProbability, MendelianAlleleProbability> getProbabilities(Ploidy father, Ploidy mother, Ploidy child) {
    if (father == Ploidy.NONE) {
      if (mother == Ploidy.NONE) {
        throw new UnsupportedOperationException("Father none, mother none, not supported.");
      }
      if (mother == Ploidy.HAPLOID) {
        if (child == Ploidy.NONE) {
          return NONE_HAPLOID_NONE;
        } else if (child == Ploidy.HAPLOID) {
          return NONE_HAPLOID_HAPLOID;
        } else {
          // Go find your real parents
          throw new UnsupportedOperationException("Mother haploid, father none, child diploid not supported.");
        }
      } else {
        throw new UnsupportedOperationException("Mother diploid, father none, not supported.");
      }
    } else if (mother == Ploidy.NONE) {
      if (father == Ploidy.HAPLOID) {
        if (child == Ploidy.NONE) {
          // haploid father, none mother, none daughter
          return HAPLOID_NONE_NONE;
        } else if (child == Ploidy.HAPLOID) {
          return HAPLOID_NONE_HAPLOID;
        } else {
          // Go find your real parents
          throw new UnsupportedOperationException("Father haploid, mother none, child diploid not supported.");
        }
      } else {
        throw new UnsupportedOperationException("Father diploid, mother none, not supported.");
      }
    } else if (child == Ploidy.NONE) {
      throw new UnsupportedOperationException("Sequence must be in child if in both parents.");
    } else {
      if (father == Ploidy.HAPLOID) {
        if (mother == Ploidy.HAPLOID) {
          if (child == Ploidy.HAPLOID) {
            return HAPLOID_HAPLOID_HAPLOID;
          } else {
            return HAPLOID_HAPLOID_DIPLOID;
          }
        }
        if (child == Ploidy.HAPLOID) {
          return HAPLOID_DIPLOID_HAPLOID;
        } else {
          return HAPLOID_DIPLOID_DIPLOID;
        }
      } else if (mother == Ploidy.HAPLOID) {
        if (child == Ploidy.HAPLOID) {
          return DIPLOID_HAPLOID_HAPLOID;
        } else {
          return DIPLOID_HAPLOID_DIPLOID;
        }
      } else if (child == Ploidy.HAPLOID) {
        throw new UnsupportedOperationException("Child must be diploid if both parents are diploid.");
      } else {
        return DIPLOID_DIPLOID_DIPLOID;
      }
    }
  }

  /**
   * @param father ploidy of father
   * @param mother ploidy of mother
   * @param child ploidy of child
   * @param logRefDenovoPrior prior for de novo mutations when parents are ref
   * @param logNonRefDenovoPrior prior for de novo mutations when parents are not ref
   * @param ref ref allele code
   * @return probability table for given configuration
   */
  public abstract MendelianAlleleProbability getMendelianAlleleProbability(Ploidy father, Ploidy mother, Ploidy child, double logRefDenovoPrior, double logNonRefDenovoPrior, int ref);

  private static final class Combined extends MendelianAlleleProbabilityFactory {
    @Override
    public MendelianAlleleProbabilityCombiner getMendelianAlleleProbability(Ploidy father, Ploidy mother, Ploidy child, double logRefDenovoPrior, double logNonRefDenovoPrior, int ref) {
      final Pair<MendelianAlleleProbability, MendelianAlleleProbability> both = getProbabilities(father, mother, child);
      return new MendelianAlleleProbabilityCombiner(both.getA(), both.getB(), logRefDenovoPrior, logNonRefDenovoPrior, ref);
    }
  }

  private static final class Mendelian extends MendelianAlleleProbabilityFactory {
    @Override
    public MendelianAlleleProbability getMendelianAlleleProbability(Ploidy father, Ploidy mother, Ploidy child, double logRefDenovoPrior, double logNonRefDenovoPrior, int ref) {
      return getProbabilities(father, mother, child).getA();
    }
  }

  private static final class Denovo extends MendelianAlleleProbabilityFactory {

    @Override
    public MendelianAlleleProbability getMendelianAlleleProbability(Ploidy father, Ploidy mother, Ploidy child, double logRefDenovoPrior, double logNonRefDenovoPrior, int ref) {
      return new DenovoWrapper(logRefDenovoPrior, logNonRefDenovoPrior, ref, getProbabilities(father, mother, child).getB());
    }
  }

  private static final class DenovoWrapper extends MendelianAlleleProbabilityDeNovo {
    private final double mDenovoLogRefPrior;
    private final double mDenovoLogNonRefPrior;
    private final int mRef;
    private final MendelianAlleleProbability mProb;

    private DenovoWrapper(double denovoLogRefPrior, double denovoLogNonRefPrior, int ref, MendelianAlleleProbability prob) {
      mDenovoLogRefPrior = denovoLogRefPrior;
      mDenovoLogNonRefPrior = denovoLogNonRefPrior;
      mRef = ref;
      mProb = prob;
    }

    @Override
    public double probabilityLn(Code code, int father, int mother, int child) {
      if (father == mRef && mother == mRef) {
        return mProb.probabilityLn(code, father, mother, child) + mDenovoLogRefPrior;
      }
      return mProb.probabilityLn(code, father, mother, child) + mDenovoLogNonRefPrior;
    }
  }
}
