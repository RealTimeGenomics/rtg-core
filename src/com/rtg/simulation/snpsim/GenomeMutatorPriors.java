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
package com.rtg.simulation.snpsim;

import java.util.Arrays;

import com.rtg.reference.Ploidy;
import com.rtg.simulation.SimulationUtils;
import com.rtg.simulation.snpsim.Mutation.DifferentMode;
import com.rtg.simulation.snpsim.Mutation.MutationType;
import com.rtg.util.PortableRandom;
import com.rtg.util.StringUtils;
import com.rtg.variant.GenomePriorParams;

/**
 * Reads and sets thresholds for genome mutator
 */
public class GenomeMutatorPriors {

  private final GenomePriorParams mPriors;
  private final Ploidy mPloid;
  /** rate for all mutations */
  private final double mRate;

  private final double mSnpThres;
  private final double mMnpThres;
  private final double mInsertThres;
  private final double mDeleteThres;
  private final double mInsDelThres;

  private final double[] mSnpModeThres;
  private final double[] mMnpModeThres;
  private final double[] mInsertModeThres;
  private final double[] mDeleteModeThres;

  private static final double VARIANT_DIFFERENT_FACTOR = 0.33333333333;

  private final double[] mMnpOmoDist;
  private final double[] mMnpEteroDist;
  private final double[] mIndelOmoDist;
  private final double[] mIndelEteroDist;
  /**
   * Array of SNP probability threshold arrays for each reference (A C G T)
   * Each array will consist of the accumulated probability of each possible call.  Calls are in the following order:
   * A C G T A:C A:G A:T C:G C:T G:T
   */
  private final double[][] mSnpThresholds = new double[4][];

  // for setting the rates from the command line
  static final double[] DEFAULT_INDEL_DIST = {0.00, 0.48, 0.17, 0.05, 0.09, 0.02, 0.03, 0.01, 0.02, 0.01, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.005, 0.015};

//  private int mStatsNumTotal;
//  private int mStatsNumEteroTotal;
//  private int[] mStatsNum;
//  private int[] mStatsNumEtero;
//
//  private int[] mStatsSnpMnpHist;
//  private int[] mStatsIndelHist;
//
//  private int[] mStatsSnpMnpHist1;
//  private int[] mStatsIndelHist1;
//  private int[] mStatsSnpMnpHist2;
//  private int[] mStatsIndelHist2;

  //                                                        SNP     ins  delete NOP
  static final double[] THRESHOLD_ENDS      = {0.5,  0.75,  1,     1};
  static final double[] THRESHOLD_ENDS_SNPS = {1,    1,     1,     1};
  static final double[] THRESHOLD_MID1      = {0.05, 0.1,   0.15,  1};
  static final double[] THRESHOLD_MID2      = {0.5,  0.7,   0.9,   1};

  /**
   * Constructor for mutation rates and threshold from priorities file
   * @param priors with priors probability
   * @param ploidy how many copies of each chromosome should be generated
   */
  public GenomeMutatorPriors(final GenomePriorParams priors, Ploidy ploidy) {
    assert priors != null;
    mPriors = priors;
    mPloid = ploidy;
    double rate = 0;
    final double[] typeRates = new double[5];
    double[] hRates = new double[3];
    double[] t;

    hRates[0] = mPriors.genomeSnpRate(false);
    hRates[1] = mPriors.genomeSnpRate(true) * (1 - VARIANT_DIFFERENT_FACTOR);
    hRates[2] = mPriors.genomeSnpRate(true) * VARIANT_DIFFERENT_FACTOR;
    mSnpModeThres = setThresholds(hRates);
    typeRates[0] = hRates[0] + hRates[1] + hRates[2]; // 0 = snp
    rate += typeRates[0];

    hRates[0] = mPriors.genomeMnpEventRate(false);
    hRates[1] = mPriors.genomeMnpEventRate(true) * (1 - VARIANT_DIFFERENT_FACTOR);
    hRates[2] = mPriors.genomeMnpEventRate(true) * VARIANT_DIFFERENT_FACTOR;
    mMnpModeThres = setThresholds(hRates);
    typeRates[1] = hRates[0] + hRates[1] + hRates[2]; // 1 = mnp
    rate += typeRates[1];

    // now it gets a bit complicated; we have to divide and conquer
    hRates = getIndelRates(mPriors);
    mInsertModeThres = setThresholds(hRates);
    typeRates[2] = hRates[0] + hRates[1] + hRates[2]; // 2 = insert
    rate += typeRates[2];

    mDeleteModeThres = mInsertModeThres;
    typeRates[3] = typeRates[2]; // 3 = delete
    rate += typeRates[3];

    if (mPloid != Ploidy.HAPLOID) {
      typeRates[4] = getInsDelRate(mPriors); // 4 == insdel
    } else {

      typeRates[4] = 0;
    }
    // has no thresholds
    rate += typeRates[4];

    mRate = rate;
    t = setThresholds(typeRates);
    mSnpThres = t[0];
    mMnpThres = t[1];
    mInsertThres = t[2];
    mDeleteThres = t[3];
    mInsDelThres = t[4];
    mMnpOmoDist = accumDistribution(mPriors.genomeMnpDistribution());
    mMnpEteroDist = accumDistribution(mPriors.genomeMnpDistribution());
    mIndelOmoDist = accumDistribution(mPriors.genomeIndelDistribution());
    mIndelEteroDist = accumDistribution(mPriors.genomeIndelDistribution(), 1);
    final double[][] snpProbs = snpProbabilities(); //Will use numbers from mPriors eventually
    for (int i = 0; i < mSnpThresholds.length; i++) {
      mSnpThresholds[i] = accumDistribution(snpProbs[i]);
    }
//    initializeStatistics();
  }

  /**
   * Generates only homozygous mutations.
   * @param snpR rate of SNPs for all mutations (0.0 - 1.0)
   * @param mnpR rate of MNPs for all mutations (0.0 - 1.0)
   * @param indelR rate of indels for all mutations (0.0 - 1.0)
   * @param totalRate total rate of mutations (per source bp)
   */
  public GenomeMutatorPriors(final Double snpR, final Double mnpR, final Double indelR, final double totalRate) {
    mPloid = Ploidy.HAPLOID;
    mPriors = null;
    double rate = 0;
    final double[] typeRates = new double[5];
    double[] t;
    double snpRate = 0;
    double mnpRate = 0;
    double indelRate = 0;

    // fill missing rates
    double remaining = totalRate;
    int split = 3;
    if (snpR != null) {
      remaining -= snpR;
      split--;
      snpRate = snpR;
    }
    if (mnpR != null) {
      remaining -= mnpR;
      split--;
      mnpRate = mnpR;
    }
    if (indelR != null) {
      remaining -= indelR;
      split--;
      indelRate = indelR;
    }

    final boolean useDefaults = snpRate == 0 && mnpRate == 0 && indelRate == 0;
    final double defValue = useDefaults ? remaining / 3 : remaining / split;
    if (useDefaults || snpR == null) {
      snpRate = defValue;
    }
    if (useDefaults || mnpR == null) {
      mnpRate = defValue;
    }
    if (useDefaults || indelR == null) {
      indelRate = defValue;
    }

    // define overall rate and thresholds
    // only homozygous simulation possible
    mSnpModeThres = new double[]{0.33, 0.66, 1.0};
    typeRates[0] = snpRate; // 0 = snp
    rate += typeRates[0];

    mMnpModeThres = new double[]{0.33, 0.66, 1.0};
    typeRates[1] = mnpRate; // 1 = mnp
    rate += typeRates[1];

    mInsertModeThres = new double[]{0.33, 0.66, 1.0};
    typeRates[2] = 0.50 * indelRate; // 2 = insert
    rate += typeRates[2];

    mDeleteModeThres = new double[]{0.33, 0.66, 1.0};
    typeRates[3] = 0.50 * indelRate; // 3 = delete
    rate += typeRates[3];

    typeRates[4] = 0.0; // 4 = insdel
    t = setThresholds(typeRates);
    rate += typeRates[4];

    mRate = rate;
    mSnpThres = t[0];
    mMnpThres = t[1];
    mInsertThres = t[2];
    mDeleteThres = t[3];
    mInsDelThres = t[4];

    mMnpOmoDist = accumDistribution(DEFAULT_INDEL_DIST);
    mMnpEteroDist = accumDistribution(DEFAULT_INDEL_DIST);
    mIndelOmoDist = accumDistribution(DEFAULT_INDEL_DIST);
    mIndelEteroDist = accumDistribution(DEFAULT_INDEL_DIST);

//    initializeStatistics();
  }

  protected double[][] snpProbabilities() {
    final double[][] probabilities = new double[4][];
    for (int i = 0; i < probabilities.length; i++) {
      probabilities[i] = new double[9];
    }
    final String[] calls = {"A", "C", "G", "T", "A:C", "A:G", "A:T", "C:G", "C:T", "G:T"};
    int a = 0, c = 0, g = 0, t = 0;
    for (int i = 0; i < calls.length; i++, a++, c++, g++, t++) {
      final double[] probs = mPriors.getPriorDistr(calls[i]);
      if (i != 0) {
        probabilities[0][a] = probs[0];
      } else {
        a--;
      }
      if (i != 1) {
        probabilities[1][c] = probs[1];
      } else {
        c--;
      }
      if (i != 2) {
        probabilities[2][g] = probs[2];
      } else {
        g--;
      }
      if (i != 3) {
        probabilities[3][t] = probs[3];
      } else {
        t--;
      }
    }
    return probabilities;
  }

  /**
   * transform prior rates into simulation indel rates
   * for delete and insert
   * @param priors the prior given
   * @return rates ordered homozygous, heterozygous with one same and different
   */
  public static double[] getIndelRates(final GenomePriorParams priors) {
    final double[] rates = new double[3];
    final double indelOmoRate = priors.genomeIndelEventRate(false);
    // homozygous rate
    rates[0] = indelOmoRate / 2.0;
    final double indelEteroRate = priors.genomeIndelEventRate(true);
    final double indelDifferentRate = indelEteroRate * VARIANT_DIFFERENT_FACTOR;
    // different rate;
    rates[2] = indelDifferentRate / 4.0;
    final double indelOneSameRate = indelEteroRate - indelDifferentRate;
    // one same rate
    rates[1] = indelOneSameRate / 2.0;
    return rates;
  }

  /**
   * Returns the rate of indels with insertion and deletion mixed
   * which are always heterozygous and different
   * @param priors the prior given
   * @return insertion-deletion rate
   */
  public static double getInsDelRate(final GenomePriorParams priors) {
    final double indelDifferentRate = priors.genomeIndelEventRate(true)
    * VARIANT_DIFFERENT_FACTOR;
    return indelDifferentRate  / 2.0;
  }

  /** @return rate of total mutations */
  public double rate() {
    return mRate;
  }

  /**
   * Used for testing
   * @param type mutation type
   * @param heterozygous if heterozygous
   * @return maximum length
   */
  public int maxLength(final MutationType type, final boolean heterozygous) {
    int len;
    if (heterozygous) {
      if (type == MutationType.MNP) {
        len = mMnpEteroDist.length;
      } else if (type == MutationType.INSERT || type == MutationType.DELETE || type == MutationType.INSDEL) {
        len = mIndelEteroDist.length + 1;
      } else {
        throw new IllegalStateException("Unpossible");
      }
    } else {
      if (type == MutationType.MNP) {
        len = mMnpOmoDist.length;
      } else if (type == MutationType.INSERT || type == MutationType.DELETE || type == MutationType.INSDEL) {
        len = mIndelOmoDist.length + 1;
      } else {
        throw new IllegalStateException("Unpossible");
      }

    }
    return len - 1;
  }

  /**
   * Sets the thresholds using an array of rates
   * @param rates rates to compute thresholds from
   * @return thresholds computed from distribution
   */
  public static double[] setThresholds(final double... rates) {
    double sum = 0;
    for (final double r : rates) {
      sum += r;
    }
    final double factor = 1.0 / sum;

    // add rates to make thresholds and normalize to 1.0
    final double[] thres = new double[rates.length];
    double currentThres = 0;
    for (int i = 0; i < rates.length; i++) {
      currentThres += rates[i];
      thres[i] = currentThres * factor;
    }
    // last one should be close to 1.0
    return thres;
  }

  private double[] accumDistribution(final double[] dist) {
    final double[] accumDist = new double[dist.length];
    double accum = 0.0;
    for (int i = 0; i < dist.length; i++) {
      accum += dist[i];
      accumDist[i] = accum;
    }
    return accumDist;
  }

  /**
   * Generates an accumulated distribution from a probability distribution
   * @param dist the non-accumulated distribution
   * @param start start point if the distribution should start at 0 or higher
   * @return the accumulated distribution
   */
  public static double[] accumDistribution(final double[] dist, final int start) {
    final double[] accumDist = new double[dist.length];
    double accum = 0.0;
    for (int i = start; i < dist.length; i++) {
      accum += dist[i];
      accumDist[i] = accum;
    }
    for (int i = start; i < dist.length; i++) {
      accumDist[i] /= accum;
    }
    return accumDist;
  }

  /**
   * Choose a length for a mutation type and depending on heterozygous
   * @param r random generator
   * @param type mutation type
   * @param heterozygous if heterozygous
   * @return next random length
   */
  public int chooseLength(final PortableRandom r, final MutationType type, final boolean heterozygous) {
    final double rand = r.nextDouble();
    return chooseLength(rand, type, heterozygous);
  }

  /**
   * Choose a length for a mutation type and depending on heterozygous
   * Used for testing
   * @param rand random value
   * @param type mutation type
   * @param heterozygous if heterozygous
   * @return length using the distribution for the given type
   */
  public int chooseLength(final double rand, final MutationType type, final boolean heterozygous) {
    if (heterozygous) {
      if (type == MutationType.MNP) {
        return SimulationUtils.chooseLength(mMnpEteroDist, rand);
      } else if (type == MutationType.INSERT || type == MutationType.DELETE || type == MutationType.INSDEL) {
        return SimulationUtils.chooseLength(mIndelEteroDist, rand) + 1;
      } else {
        throw new IllegalStateException("Unpossible");
      }
    } else {
      if (type == MutationType.MNP) {
        return SimulationUtils.chooseLength(mMnpOmoDist, rand);
      } else if (type == MutationType.INSERT || type == MutationType.DELETE || type == MutationType.INSDEL) {
        return SimulationUtils.chooseLength(mIndelOmoDist, rand) + 1; // the indel distribution starts from length 1
      } else {
        throw new IllegalStateException("Unpossible");
      }
    }
  }

  /**
   * Choose a mutation type
   * @param r random generator
   * @return mutation type
   */
  public MutationType chooseType(final PortableRandom r) {
    final double rand = r.nextDouble();
    return chooseType(rand);
  }

  /**
   * Choose a mutation type
   * @param rand value generated for selection
   * @return mutation type
   */
  public MutationType chooseType(final double rand) {
    if (rand < mSnpThres) {
      return MutationType.SNP;
    } else if (rand < mMnpThres) {
      return MutationType.MNP;
    } else if (rand < mInsertThres) {
      return MutationType.INSERT;
    } else if (rand < mDeleteThres) {
      return MutationType.DELETE;
    } else if (rand < mInsDelThres) {
      return MutationType.INSDEL;
    }
    throw new IllegalArgumentException("Invalid mutation distribution");
  }

  /**
   * Chooses if homozygous, heterozygous with one same or
   * or heterozygous with both different
   * @param r random generator
   * @param type type of mutation
   * @return mode of mutation
   */
  public DifferentMode chooseDifferentMode(final PortableRandom r, final MutationType type) {
    final double rand = r.nextDouble();
    return chooseDifferentMode(rand, type);
  }

  /**
   * Chooses mode of mutation, interface used for testing
   * @param rand random double
   * @param type type of mutation
   * @return mode of mutation
   */
  public DifferentMode chooseDifferentMode(final double rand, final MutationType type) {
    if (mPloid == Ploidy.HAPLOID) {
      return DifferentMode.HOMOZYGOUS;
    }
    switch (type.simple()) {
    case MNP:
      return chooseDifferentMode(rand, mMnpModeThres);
    case SNP:
      return chooseDifferentMode(rand, mSnpModeThres);
    case INSERT:
      return chooseDifferentMode(rand, mInsertModeThres);
    case DELETE:
      return chooseDifferentMode(rand, mDeleteModeThres);
    case INSDEL:
      return DifferentMode.DIFFERENT;
    default:
      throw new IllegalStateException("Unpossible");
    }
  }

  private DifferentMode chooseDifferentMode(final double rand, final double[] thres) {
    if (rand < thres[0]) {
      return DifferentMode.HOMOZYGOUS;
    } else if (rand < thres[1]) {
      return DifferentMode.ONE_ONLY;
    }
    return DifferentMode.DIFFERENT;
  }

  protected boolean hasPriors() {
    return mPriors != null;
  }

  /**
   * Return a detailed human readable representation of this object.
   * It is intended that this show internal details of the
   * object structure that may be relevant to an implementor/debugger but not to
   * a user.
   *
   * @return string with object representation.
   */
  @Override
  public String toString() {
    return ("Rate " + mRate + StringUtils.LS + "SnpThres " + mSnpThres + StringUtils.LS) + "MnpThres " + mMnpThres + StringUtils.LS + "InsertThres " + mInsertThres + StringUtils.LS + "DeleteThres " + mDeleteThres + StringUtils.LS + "InsDelThres " + mInsDelThres + StringUtils.LS + "SnpEteroThres " + mSnpModeThres[0] + StringUtils.LS + "MnpEteroThres " + mMnpModeThres[0] + StringUtils.LS + "InsertEteroThres " + mInsertModeThres[0] + StringUtils.LS + "DeleteEteroThres " + mDeleteModeThres[0] + StringUtils.LS + "InsDelEteroThres n/a" + StringUtils.LS + "VariantDifferentFactor " + VARIANT_DIFFERENT_FACTOR + StringUtils.LS + "MnpOmoDist length " + mMnpOmoDist.length + StringUtils.LS + "MnpOmoDist: " + Arrays.toString(mMnpOmoDist) + StringUtils.LS + "MnpEteroDist length " + mMnpEteroDist.length + StringUtils.LS + "MnpEteroDist: " + Arrays.toString(mMnpEteroDist) + StringUtils.LS + "IndelOmoDist length " + mIndelOmoDist.length + StringUtils.LS + "IndelOmoDist: " + Arrays.toString(mIndelOmoDist) + StringUtils.LS + "IndelEteroDist length " + mIndelEteroDist.length + StringUtils.LS + "IndelEteroDist: " + Arrays.toString(mIndelEteroDist) + StringUtils.LS + "MnpThresholds-ends: " + Arrays.toString(THRESHOLD_ENDS) + StringUtils.LS + "MnpThresholds-ends-snponly: " + Arrays.toString(THRESHOLD_ENDS_SNPS) + StringUtils.LS + "MnpThresholds-mid1: " + Arrays.toString(THRESHOLD_MID1) + StringUtils.LS + "MnpThresholds-mid2: " + Arrays.toString(THRESHOLD_MID2) + StringUtils.LS;
  }

  /**
   * @param r random generator
   * @param refBase the ordinal of the reference nucleotide
   * @return the pair of generated nucleotides
   */
  public byte[] chooseAltSnp(PortableRandom r, byte refBase) {
    return chooseAltSnp(r.nextDouble(), refBase);
  }

  /**
   * @param rand random double
   * @param refBase the ordinal of the reference nucleotide
   * @return the pair of generated nucleotides
   */
  public byte[] chooseAltSnp(double rand, byte refBase) {
    if (refBase == 0) {
      //Special case, have asked for SNP on N nucleotide
      return new byte[] {0, 0};
    }
    //A C T G A:C A:T A:G C:T C:G G:T
    final double[] thresholds = mSnpThresholds[refBase - 1];
    final double adjustedRand = rand * thresholds[8];
    if (adjustedRand < thresholds[0]) {
      if (refBase == 1) {
        return new byte[] {2, 2};
      } else {
        return new byte[] {1, 1};
      }
    } else if (adjustedRand < thresholds[1]) {
      if (refBase <= 2) {
        return new byte[] {3, 3};
      } else {
        return new byte[] {2, 2};
      }
    } else if (adjustedRand < thresholds[2]) {
      if (refBase <= 3) {
        return new byte[] {4, 4};
      } else {
        return new byte[] {3, 3};
      }
    } else if (adjustedRand < thresholds[3]) {
      return new byte[] {1, 2};
    } else if (adjustedRand < thresholds[4]) {
      return new byte[] {1, 3};
    } else if (adjustedRand < thresholds[5]) {
      return new byte[] {1, 4};
    } else if (adjustedRand < thresholds[6]) {
      return new byte[] {2, 3};
    } else if (adjustedRand < thresholds[7]) {
      return new byte[] {2, 4};
    } else if (adjustedRand < thresholds[8]) {
      return new byte[] {3, 4};
    }
    throw new IllegalArgumentException("Invalid snp distribution");
  }
}
