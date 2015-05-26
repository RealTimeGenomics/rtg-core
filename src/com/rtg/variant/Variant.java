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

import java.util.HashSet;
import java.util.Set;

import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.NoNonIdentityMeasure;

/**
 * Container for Variance Locus, parent Call, and n sample Calls.
 * Singleton implemented as a single sample.
 * Parent call used to generate Info field in VCF, and sample Calls to generate
 * sample columns.
 */
public class Variant extends IntegralAbstract implements Comparable<Variant> {

  /**
   * Enumeration of variant filter types
   */
  public enum VariantFilter {
    /** filtered due to coverage */
    COVERAGE(0x01),
    /** filtered due to ambiguity */
    AMBIGUITY(0x02),
    /** filtered due to complex equivalence */
    COMPLEX_EQUIVALENT(0x04),
    /** filtered due to hyper complex */
    HYPER_COMPLEX(0x08),
    /** filtered due to ion torrent */
    IONTORRENT(0x010),
    /** filtered due to failing complex calling */
    FAILED_COMPLEX(0x020),
    /** filtered due to falling outside bed regions */
    BED_REGION(0x040),
    /** filtered due to other */
    OTHER(0x080);

    private final int mMask;

    private VariantFilter(int mask) {
      mMask = mask;
    }

    /**
     * @return the mask of the variant filter
     */
    public int mask() {
      return mMask;
    }
  }

  private final VariantLocus mLocus;
  private int mFilters;
  private final VariantSample[] mSamples;

  private int mIndelLength;
  private boolean mIndel = false;
  private boolean mOverflow = false;
  private boolean mInvalidRef = false; // Reference did not form a hypothesis (i.e. due to N's in ref)
  private boolean mComplexScored = false;
  private boolean mComplexEquivalent = false;
  private boolean mInteresting = false;
  private boolean mTrimmed = false;
  private int mSplitId = -1;
  private boolean mForceComplex = false;

  private Double mNonIdentityPosterior = null;

  private String mPossibleCause;
  private Double mPossibleCauseScore;
  private Double mNormalCancerScore;
  private Double mDiseasePresenceScore;
  private Double mLoh;

  /**
   * @param locus the locus for this variant
   * @param filters the set of filters to apply
   * @param samples the samples for this call, null if Ploidy.NONE for that sample
   */
  public Variant(VariantLocus locus, int filters, VariantSample... samples) {
    if (locus == null) {
      throw new NullPointerException();
    }
    mLocus = locus;
    mSamples = samples;
    mFilters = filters;
    if (samples.length > 0) {
      setSampleLikelihoods(samples, locus.getRefNts());
    }
  }

  /**
   * No filters applied.
   * @param locus the locus for this variant
   * @param samples the sample for this call, null if Ploidy.NONE for that sample
   */
  public Variant(VariantLocus locus, VariantSample... samples) {
    this(locus, 0, samples);
  }

  private static void setSampleLikelihoods(VariantSample[] samples, String ref) {
    final Set<String> alleleSet = alleleSet(samples, ref);
    for (VariantSample s : samples) {
      // a sample may already have a genotype likelihood as a result of the variant splitting/trimming don't overwrite it
      if (s != null && s.getMeasure() != null && s.getGenotypeLikelihoods() == null) {
        if (!(s.getMeasure() instanceof NoNonIdentityMeasure)) {
          s.setGenotypeLikelihoods(s.computeGenotypeLikelihoods(alleleSet));
        }
      }
    }
  }

  static Set<String> alleleSet(VariantSample[] samples, String ref) {
    final Set<String> alleles = new HashSet<>();
    alleles.add(ref);
    for (final VariantSample sample : samples) {
      if (sample != null  && sample.getMeasure() != null) {
        final int best = sample.getMeasure().best();
        final Description desc = sample.getMeasure().hypotheses().description();
        final Code code = sample.getMeasure().hypotheses().code();
        alleles.add(desc.name(code.a(best)));
        alleles.add(desc.name(code.bc(best)));
      }
    }
    return alleles;
  }

  public VariantLocus getLocus() {
    return mLocus;
  }

  public int getNumberOfSamples() {
    return mSamples.length;
  }

  /**
   * Return a particular sample
   * @param index index of the sample
   * @return sample at the requested index, null if Ploidy.NONE for that sample
   */
  public VariantSample getSample(int index) {
    return mSamples[index];
  }

  public boolean isFiltered() {
    return mFilters > 0;
  }

  public int getFilters() {
    return mFilters;
  }

  /**
   * Get trimmed.
   * @return Returns the trimmed.
   */
  public boolean isTrimmed() {
    return mTrimmed;
  }

  /**
   * Set trimmed.
   */
  public void setTrimmed() {
    mTrimmed = true;
  }

  /**
   * Get split id.
   * @return Returns the split id.
   */
  public int getSplitId() {
    return mSplitId;
  }

  /**
   * Set split id.
   * @param splitId The split id to set.
   */
  public void setSplitId(int splitId) {
    mSplitId = splitId;
  }

  /**
   * @param filter a filter to check
   * @return true if this call is filtered due to the given filter
   */
  public boolean isFiltered(VariantFilter filter) {
    return (mFilters & filter.mask()) > 0;
  }

  /**
   * Add a reason why this call is filtered.
   * @param filter the filter to add
   */
  public void addFilter(VariantFilter filter) {
    if (isFiltered(filter)) {
      return;
    }
    mFilters |= filter.mask();
  }

  /** Set if this call has no valid reference hypothesis (i.e. Ref contained N's) */
  public void setInvalidRef() {
    mInvalidRef = true;
  }

  /** @return true if this call did not have a valid reference hypothesis */
  public boolean isInvalidRef() {
    return mInvalidRef;
  }

  /**
   * Set if this call is considered interesting
   */
  public void setInteresting() {
    mInteresting = true;
  }
  /**
   * @return true if this call is considered interesting
   */
  public boolean isInteresting() {
    return mInteresting;
  }

  /**
   * Set if this call has been scored by the complex caller
   */
  public void setComplexScored() {
    mComplexScored = true;
  }

  /**
   * @return true if this call has been scored by the complex caller
   */
  public boolean isComplexScored() {
    return mComplexScored;
  }

  public int getIndelLength() {
    return mIndelLength;
  }

  /**
   * Indicates that this call is an insert or delete of a specified length, to be processed by the complex caller.
   * @param length of operation in cigar
   */
  public void setIndel(int length) {
    mIndel = true;
    mForceComplex = true;
    mIndelLength = length;
  }

  /**
   * @return true if this call is an insert or delete (i.e., to be processed by the complex caller)
   */
  public boolean isIndel() {
    return mIndel;
  }

  /**
   * Mark this record as indicating a region of excessive coverage.
   */
  public void setOverflow() {
    mOverflow = true;
  }

  /**
   * @return true if this record is indicating a region of excessive coverage.
   */
  public boolean isOverflow() {
    return mOverflow;
  }

  /**
   * @return true if this record should force complex calling.
   */
  public boolean isForceComplex() {
    return mForceComplex;
  }

  /**
   * Mark this record as needing force calling to be forced.
   */
  public void setForceComplex() {
    mForceComplex = true;
  }


  /**
   * Set if this call is a complex equivalent
   */
  public void setComplexEquivalent() {
    mComplexEquivalent = true;
  }
  public boolean isComplexEquivalent() {
    return mComplexEquivalent;
  }

  public Double getNonIdentityPosterior() {
    return mNonIdentityPosterior;
  }

  public void setNonIdentityPosterior(Double posterior) {
    mNonIdentityPosterior = posterior;
  }

  /**
   * Get possible cause.
   * @return Returns the possible cause.
   */
  public String getPossibleCause() {
    return mPossibleCause;
  }

  /**
   * Set possible cause.
   * @param possibleCause The possible cause to set.
   */
  public void setPossibleCause(String possibleCause) {
    mPossibleCause = possibleCause;
  }

  /**
   * Get possible cause score.
   * @return Returns the possible cause score.
   */
  public Double getPossibleCauseScore() {
    return mPossibleCauseScore;
  }

  /**
   * Set possible cause score.
   * @param possibleCauseScore The possible cause score to set.
   */
  public void setPossibleCauseScore(Double possibleCauseScore) {
    mPossibleCauseScore = possibleCauseScore;
  }

  /**
   * Get the log posterior for the cancer call being different from the normal call.
   * @return Returns the normal cancer score.
   */
  public Double getNormalCancerScore() {
    return mNormalCancerScore;
  }

  /**
   * Set the log posterior for the cancer call being different from the normal call.
   * @param ncScore The normal cancer score to set.
   */
  public void setNormalCancerScore(Double ncScore) {
    mNormalCancerScore = ncScore;
  }

  /**
   * Get the log posterior for the presence of disease
   * @return Returns the disease presence score
   */
  public Double getDiseasePresenceScore() {
    return mDiseasePresenceScore;
  }

  /**
   * Set the log posterior for a disease presence score
   * @param diseasePresenceScore The disease presence score to set.
   */
  public void setDiseasePresenceScore(Double diseasePresenceScore) {
    mDiseasePresenceScore = diseasePresenceScore;
  }

  public void setLoh(Double loh) {
    mLoh = loh;
  }

  public Double getLoh() {
    return mLoh;
  }

  /**
   * Method to copy all non final fields from one object to another
   * should ignore copying split id
   *
   * @param copyFrom values to copy from
   * @param copyTo value to copy to
   */
  public static void copy(Variant copyFrom, Variant copyTo) {
    if (copyFrom.isComplexEquivalent()) {
      copyTo.setComplexEquivalent();
    }
    if (copyFrom.isComplexScored()) {
      copyTo.setComplexScored();
    }
    if (copyFrom.isInvalidRef()) {
      copyTo.setInvalidRef();
    }
    if (copyFrom.isIndel()) {
      copyTo.setIndel(copyFrom.mIndelLength);
    }
    if (copyFrom.isForceComplex()) {
      copyTo.setForceComplex();
    }
    if (copyFrom.isInteresting()) {
      copyTo.setInteresting();
    }
    copyTo.setPossibleCause(copyFrom.getPossibleCause());
    copyTo.setPossibleCauseScore(copyFrom.getPossibleCauseScore());
    copyTo.setNormalCancerScore(copyFrom.getNormalCancerScore());
    copyTo.setDiseasePresenceScore(copyFrom.getDiseasePresenceScore());
      copyTo.setNonIdentityPosterior(copyFrom.getNonIdentityPosterior());
    if (copyFrom.isTrimmed()) {
      copyTo.setTrimmed();
    }
    copyTo.mFilters = copyFrom.mFilters;
    copyTo.setLoh(copyFrom.getLoh());
  }

  @Override
  public int compareTo(Variant that) {
    final int diff = this.getLocus().getStart() - that.getLocus().getStart();
    if (diff != 0) {
      return diff;
    }

    final int thisO = this.isIndel() ? 0 : 1;
    final int thatO = that.isIndel() ? 0 : 1;
    return thisO - thatO;
  }

  @Override
  public boolean equals(Object arg0) {
    return arg0 instanceof Variant && compareTo((Variant) arg0) == 0;
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(getLocus().getStart(), getTypeString().hashCode());
  }

  /**
   * Is this call a strict SNP (i.e., not a MNP or indel).
   * @return true iff this is a SNP
   */
  public boolean isSnp() {
    if (isFiltered(VariantFilter.FAILED_COMPLEX)) {
      return false;
    }
    if (getLocus().getRefNts() == null
      || getLocus().getRefNts().length() != 1) {
      return false;
    }
    Boolean noCats = null;
    for (int i = 0; i < getNumberOfSamples(); i++) {
      if (getSample(i) == null || getSample(i).getName() == null) {
        if (noCats == null) {
          noCats = true;
        }
        continue;
      } else {
        noCats = false;
      }
      final String name = getSample(i).getName();
      final int len = name.length();
      if (len != 1 && len != 3) {
        return false;
      }
      if (len == 3 && name.charAt(1) != ':') {
        return false;
      }
    }
    if (noCats) {
      return true;
    }

    return getLocus().getEnd() - getLocus().getStart() == 1;
  }

  /** @return type of this variant. */
  String getTypeString() {
    if (isIndel()) {
      return "INDEL";
    } else {
      return "SNP";
    }
  }

  /**
   * Check that calls have names
   * @return true if the first non null sample has a name or if all samples are null
   */
  public boolean hasCallNames() {
    for (int i = 0; i < getNumberOfSamples(); i++) {
      if (getSample(i) != null) {
        return getSample(i).getName() != null; //only need to check first non-null sample, if one BC is null they all will be.
      }
    }
    return true;
  }


  @Override
  public boolean integrity() {
    Exam.assertNotNull(getLocus());
    Exam.assertTrue(getLocus().getStart() >= 0);
    Exam.assertTrue(getNumberOfSamples() > 0 || isOverflow());
    final int minLen;
    final int maxLen;
    if (isIndel()) {
      //        case END_HYPER:
      minLen = 0;
      maxLen = 0;
    } else if (isComplexScored()) {
      minLen = 0;
      maxLen = Integer.MAX_VALUE - getLocus().getStart();
    } else {
      minLen = 1;
      maxLen = 1;
    }
    Exam.assertTrue(getLocus().getEnd() >= getLocus().getStart() + minLen);
    Exam.assertTrue("end=" + getLocus().getEnd() + " start=" + getLocus().getStart() + " max=" + maxLen + " type=" + getTypeString(), getLocus().getEnd() <= getLocus().getStart() + maxLen);
    return true;
  }
}
