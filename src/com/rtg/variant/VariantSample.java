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

package com.rtg.variant;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.rtg.reference.Ploidy;
import com.rtg.variant.VariantSampleInfo.VariantFormatEnum;
import com.rtg.variant.bayes.GenotypeMeasure;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.Statistics;

/**
 * Contains details of a call, but nothing of position.
 */
public class VariantSample {

  /** enum */
  public enum DeNovoStatus {
    /** call is detected as de novo on this sample */
    IS_DE_NOVO,
    /** call is detected as not de novo on this sample*/
    NOT_DE_NOVO,
    /** not enough information is available to determine if call is de novo, or no de novo calls are found at site */
    UNSPECIFIED
  }

  //immutable
  private final String mName;
  private final boolean mIsIdentity;
  private final Ploidy mPloidy;

  private Statistics<?> mStats;

  //potentially mutable
  private final VariantSampleInfo[] mInfo = new VariantSampleInfo[VariantFormatEnum.values().length];
  private final DeNovoStatus mIsDeNovo;
  private String mVariantAllele; // The variant allele, used during somatic calling

  private final GenotypeMeasure mMeasure;
  private final Double mDeNovoPosterior;

  private Map<Set<String>, Double> mGenotypeLikelihoods = null;

  /**
   * Create a variant sample
   * @param ploidy ploidy of this call
   * @param name name of hypothesis
   * @param isIdentity true if identity
   * @param measure posterior measure
   * @param isDeNovo as per semantics described on enum
   * @param deNovoPosterior posterior score that this is a de novo mutation
   */
  public VariantSample(Ploidy ploidy, String name, boolean isIdentity, GenotypeMeasure measure, DeNovoStatus isDeNovo, Double deNovoPosterior) {
    assert isDeNovo != null;
    mPloidy = ploidy;
    mName = name;
    mIsIdentity = isIdentity;
    mMeasure = measure;
    mIsDeNovo = isDeNovo;
    mDeNovoPosterior = deNovoPosterior;
    mVariantAllele = null;
  }

  /**
   * Create a variant sample with no hypothesis information. Generally used for error conditions i.e. over coverage
   * @param ploidy ploidy of this call
   */
  public VariantSample(Ploidy ploidy) {
    this(ploidy, null, true, null, DeNovoStatus.UNSPECIFIED, null);
  }

  public GenotypeMeasure getMeasure() {
    return mMeasure;
  }

  public Ploidy getPloidy() {
    return mPloidy;
  }

  /**
   * @return the name of the genotype hypothesis
   */
  public String getName() {
    return mName;
  }

  public boolean isIdentity() {
    return mIsIdentity;
  }

  /**
   * @return the posterior score for the hypothesis (normalized in natural log space and of form <code>ln(p/(1-p))</code>).
   */
  public Double getPosterior() {
    if (mMeasure == null) {
      return 0.0; // ln(1)
    }
    return mMeasure.arithmetic().poss2Ln(mMeasure.bestPosterior());
  }

  /**
   * @return the non identity posterior posterior or null if there is none
   */
  public Double getNonIdentityPosterior() {
    // TODO Disease caller doesn't produce a non identity posterior score so check for NaN
    if (mMeasure == null || Double.isNaN(mMeasure.nonIdentityPosterior())) {
      return null;
    }
    return mMeasure.arithmetic().poss2Ln(mMeasure.nonIdentityPosterior());
  }

  public DeNovoStatus isDeNovo() {
    return mIsDeNovo;
  }

  public Double getDeNovoPosterior() {
    return mDeNovoPosterior;
  }

  public String getVariantAllele() {
    return mVariantAllele;
  }

  public void setVariantAllele(String name) {
    mVariantAllele = name;
  }

  public Statistics<?> getStats() {
    return mStats;
  }

  public void setStats(Statistics<?> stats) {
    mStats = stats;
  }

  /**
   * Set the coverage and coverage correction values
   * @param coverage the coverage of this call
   * @param coverageCorrection the coverage correction value for this call
   */
  public void setCoverage(int coverage, double coverageCorrection) {
    setCoverage(coverage);
    setCoverageCorrection(coverageCorrection);
  }

  /**
   * Set the coverage value
   * @param coverage the coverage of this call
   */
  public void setCoverage(int coverage) {
    if (mInfo[VariantFormatEnum.COVERAGE.ordinal()] == null) {
      mInfo[VariantFormatEnum.COVERAGE.ordinal()] = new VariantSampleInfo(coverage);
    } else {
      throw new RuntimeException("Only call this once or it won't work");
    }
  }

  /**
   * Set the coverage correction value
   * @param coverageCorrection the coverage correction value for this call
   */
  public void setCoverageCorrection(double coverageCorrection) {
    setDoubleInfoHelper(VariantFormatEnum.COVERAGE_CORRECTION, coverageCorrection);
  }

  public Integer getCoverage() {
    return mInfo[VariantFormatEnum.COVERAGE.ordinal()] == null ? null : mInfo[VariantFormatEnum.COVERAGE.ordinal()].intValue();
  }
  public Double getCorrection() {
    return getDoubleInfoHelper(VariantFormatEnum.COVERAGE_CORRECTION);
  }

 /**
  * Get ambiguity ratio.
  * @return Returns the ambiguity ratio.
  */
 public Double getAmbiguityRatio() {
   return getDoubleInfoHelper(VariantFormatEnum.AMBIGUITY_RATIO);
 }

 /**
  * Set ambiguity ratio.
  * @param ambiguityRatio The ambiguity ratio to set.
  */
 public void setAmbiguityRatio(Double ambiguityRatio) {
   setDoubleInfoHelper(VariantFormatEnum.AMBIGUITY_RATIO, ambiguityRatio);
 }

  /**
   * @param val strand bias for allele one
   */
  public void setHoeffdingStrandBiasAllele1(Double val) {
    setDoubleInfoHelper(VariantFormatEnum.HOEFFDING_STRAND_BIAS_A1, val);
  }

  /**
   * @param val strand bias for allele two
   */
  public void setHoeffdingStrandBiasAllele2(Double val) {
    setDoubleInfoHelper(VariantFormatEnum.HOEFFDING_STRAND_BIAS_A2, val);
  }

  /**
   * @param val allele balance for heterozygous calls
   */
  public void setHoeffdingAlleleBalanceHet(Double val) {
    setDoubleInfoHelper(VariantFormatEnum.HOEFFDING_ALLELE_BALANCE_HET, val);
  }

  /**
   * @param val unmated bias for allele one
   */
  public void setHoeffdingUnmatedBiasAllele1(Double val) {
    setDoubleInfoHelper(VariantFormatEnum.HOEFFDING_UNMATED_BIAS_A1, val);
  }

  /**
   * @param val unmated bias for allele two
   */
  public void setHoeffdingUnmatedBiasAllele2(Double val) {
    setDoubleInfoHelper(VariantFormatEnum.HOEFFDING_UNMATED_BIAS_A2, val);
  }

  /**
   * @param val allele balance for homozygous calls
   */
  public void setHoeffdingAlleleBalanceHom(Double val) {
    setDoubleInfoHelper(VariantFormatEnum.HOEFFDING_ALLELE_BALANCE_HOM, val);
  }

  public Double getHoeffdingStrandBiasAllele1() {
    return getDoubleInfoHelper(VariantFormatEnum.HOEFFDING_STRAND_BIAS_A1);
  }

  public Double getHoeffdingStrandBiasAllele2() {
    return getDoubleInfoHelper(VariantFormatEnum.HOEFFDING_STRAND_BIAS_A2);
  }

  public Double getHoeffdingAlleleBalanceHet() {
    return getDoubleInfoHelper(VariantFormatEnum.HOEFFDING_ALLELE_BALANCE_HET);
  }
  public Double getHoeffdingAlleleBalanceHom() {
    return getDoubleInfoHelper(VariantFormatEnum.HOEFFDING_ALLELE_BALANCE_HOM);
  }

  public Double getHoeffdingUnmatedBiasAllele1() {
    return getDoubleInfoHelper(VariantFormatEnum.HOEFFDING_UNMATED_BIAS_A1);
  }

  public Double getHoeffdingUnmatedBiasAllele2() {
    return getDoubleInfoHelper(VariantFormatEnum.HOEFFDING_UNMATED_BIAS_A2);
  }

  /**
   * @param val placed unmapped ratio
   */
  public void setPlacedUnmappedRatio(Double val) {
    setDoubleInfoHelper(VariantFormatEnum.UNPLACED_RATIO, val);
  }

  public Double getPlacedUnmappedRatio() {
    return getDoubleInfoHelper(VariantFormatEnum.UNPLACED_RATIO);
  }

  /**
   * @param val read position correct probability
   */
  public void setHoeffdingReadPositionBias(Double val) {
    setDoubleInfoHelper(VariantFormatEnum.HOEFFDING_READ_POSITION, val);
  }

  public Double getHoeffdingReadPositionBias() {
    return getDoubleInfoHelper(VariantFormatEnum.HOEFFDING_READ_POSITION);
  }

  public String getStatisticsString() {
    return mInfo[VariantFormatEnum.STATISTICS.ordinal()] == null ? null : mInfo[VariantFormatEnum.STATISTICS.ordinal()].stringValue();
  }
  public void setStatisticsString(String statisticsString) {
    mInfo[VariantFormatEnum.STATISTICS.ordinal()] = new VariantSampleInfo(statisticsString);

  }
  /**
   * Append additional statistics to the statistics string
   * @param string the string to add
   */
  public void appendStatisticsString(String string) {
    if (mInfo[VariantFormatEnum.STATISTICS.ordinal()] == null) {
      setStatisticsString(string);
    } else {
      setStatisticsString(mInfo[VariantFormatEnum.STATISTICS.ordinal()].stringValue() + string);
    }
  }

  private Double getDoubleInfoHelper(VariantFormatEnum info) {
    return mInfo[info.ordinal()] == null ? null : mInfo[info.ordinal()].doubleValue();
  }

  private void setDoubleInfoHelper(VariantFormatEnum info, Double value) {
    assert mInfo[info.ordinal()] == null;
    if (value != null) {
      mInfo[info.ordinal()] = new VariantSampleInfo(value);
    }
  }

  private Map<Set<String>, Double> computeGenotypeLikelihoodsHaploid(Set<String> alleles) {
    final Map<Set<String>, Double> result = new HashMap<>(alleles.size());
    final Hypotheses<?> h = mMeasure.hypotheses();
    for (String allele : alleles) {
      final int index = h.description().indexOf(allele);
      final double poss;
      if (index == -1) {
        return Collections.emptyMap();
      } else {
        poss = mMeasure.measure(index);
      }
      result.put(Collections.singleton(allele), mMeasure.arithmetic().poss2Ln(poss));
    }
    return result;
  }

  private Map<Set<String>, Double> computeGenotypeLikelihoodsDiploid(Set<String> alleles) {
    final Map<Set<String>, Double> result = new HashMap<>();
    final List<String> alleleList = new ArrayList<>(alleles);
    final Hypotheses<?> h = mMeasure.hypotheses();
    for (int i = 0; i < alleleList.size(); ++i) {
      final String a = alleleList.get(i);
      final int aIndex = h.description().indexOf(a);
      for (int j = i; j < alleleList.size(); ++j) {
        final String b = alleleList.get(j);
        final int bIndex = h.description().indexOf(b);
        final double poss;
        if (aIndex == -1 || bIndex == -1) {
          return Collections.emptyMap();
        } else {
          final int hypIndex = h.code().code(aIndex, bIndex);
          poss = mMeasure.measure(hypIndex);
        }
        result.put(pairSet(a, b), mMeasure.arithmetic().poss2Ln(poss));
      }
    }
    return result;
  }

  /**
   * Produce a set containing the pair of strings
   * @param a first string
   * @param b second string
   * @return a set containing a &amp; b
   */
  public static Set<String> pairSet(String a, String b) {
    final Set<String> set = new HashSet<>();
    set.add(a);
    set.add(b);
    return set;
  }

  /**
   * Given a set of alleles compute what your likelihoods should be
   * @param alleles set of alleles to produce likelihoods for
   * @return Map from allele combination to likelihood
   */
  public Map<Set<String>, Double> computeGenotypeLikelihoods(Set <String> alleles) {
    final Hypotheses<?> hyp = mMeasure.hypotheses();
    if (alleles.size() < 1)  {
      return Collections.emptyMap();
    } else if (hyp.ploidy() == Ploidy.HAPLOID) {
      return computeGenotypeLikelihoodsHaploid(alleles);
    } else if (hyp.ploidy() == Ploidy.DIPLOID) {
      return computeGenotypeLikelihoodsDiploid(alleles);
    } else if (hyp.ploidy() == Ploidy.NONE) {
      return Collections.emptyMap();
    } else {
      throw new UnsupportedOperationException("Cannot calculate GL field for unsupported ploidy: " + hyp.ploidy());
    }
  }

  /**
   * Accessor
   * @return a mapping from allele set to un-normalised measure
   */
  public Map<Set<String>, Double> getGenotypeLikelihoods() {
    return mGenotypeLikelihoods;
  }

  /**
   * Explicitly set the un-normalised genotype likelihoods
   * @param likelihoods map from alleles in the genotype to un-normalised likelihood
   */
  public void setGenotypeLikelihoods(Map<Set<String>, Double> likelihoods) {
    mGenotypeLikelihoods = likelihoods;
  }


  @Override
  public String toString() {
    return getName();
  }

  /**
   * Method to copy all non final fields from one object to another
   * @param copyFrom values to copy from
   * @param copyTo value to copy to
   */
  public static void copy(VariantSample copyFrom, VariantSample copyTo) {
    for (int i = 0; i < copyFrom.mInfo.length; ++i) {
      if (copyFrom.mInfo[i] != null) {
        copyTo.mInfo[i] = copyFrom.mInfo[i];
      }
    }
    copyTo.setVariantAllele(copyFrom.getVariantAllele());
    copyTo.setGenotypeLikelihoods(copyFrom.getGenotypeLikelihoods());
    copyTo.setStats((Statistics<?>) copyFrom.getStats().copy());
  }
}
