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

package com.rtg.variant.bayes.multisample.cancer;

import java.util.Arrays;

import com.rtg.relation.Relationship;
import com.rtg.relation.Relationship.RelationshipType;
import com.rtg.util.StringUtils;
import com.rtg.util.format.FormatReal;
import com.rtg.variant.PerSampleVariantStatistics;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantStatistics;
import com.rtg.variant.format.VcfFormatField;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.header.VcfHeader;

/**
 * Specialization of variant statistics for the cancer caller.  It only outputs our usual ratio
 * statistics for the derived sample.
 *
 * Tracking and computation of numbers needed to make an estimate of contamination in the cancer sample.
 * We are keeping track of statistics for two samples.  The normal sample is tracked so that we can
 * estimate the reference bias (which is assumed to carry over to the cancer sample).  For this
 * we want to keep the best looking 0/1 normal variants.  The cancer samples is tracked so that we
 * can estimate the contamination of the cancer sample with normal.  This is done by recording the
 * best 0/0 -&gt; 0/1 cancer calls.  The contamination estimate also takes the normal reference bias
 * into account.  In both cases, the same ranking field is used which is assumed to be of type double
 * (with null for missing values) and where bigger is better.
 */
public class SomaticStatistics extends VariantStatistics {

  private static final String GT = VcfFormatField.GT.name();
  private static final String AD = VcfFormatField.AD.name();
  private static final int[] DIPLOID00 = {0, 0};
  private static final int[] DIPLOID01 = {0, 1};
  private static final int BASIS = 1000;
  private static final FormatReal FORMAT_REAL = new FormatReal(3, 2);

  private final String mRankingField;
  private final String mNormalSampleName;
  private final String mCancerSampleName;
  private final TopScorer mNormalStore;
  private final TopScorer mCancerStore;

  /**
   * Somatic specific variant statistics which only outputs the statistics for the derived sample.
   * @param params the params
   * @param rankingField field to use for ranking variants
   */
  public SomaticStatistics(VariantParams params, String rankingField) {
    super(params.directory());
    final Relationship[] derived = params.genomeRelationships().relationships(RelationshipType.ORIGINAL_DERIVED);
    final int numberOfGenomes = params.genomeRelationships().genomes().length;
    assert derived.length == 1 || numberOfGenomes == 2;
    final String[] genomeNames = {derived[0].first(), derived[0].second()};
    mRankingField = rankingField;
    mNormalSampleName = genomeNames[0];
    mCancerSampleName = genomeNames[1];
    mNormalStore = new TopScorer(BASIS);
    mCancerStore = new TopScorer(BASIS);
  }

  private boolean isGenotype(final VcfRecord rec, final int sampleId, final int[] desired) {
    // desired is assumed to be sorted
    final int[] gt = VcfUtils.splitGt(rec.getSampleString(sampleId, GT));
    if (gt.length == desired.length) {
      Arrays.sort(gt);
      return Arrays.equals(gt, desired);
    }
    return false;
  }

  private void insert(final TopScorer store, final Double rankingScore, final String alleleDepth) {
    final String[] alleleCounts = StringUtils.split(alleleDepth, ',');
    store.add(rankingScore, Integer.parseInt(alleleCounts[0]), Integer.parseInt(alleleCounts[1]));
  }

  @Override
  public void tallyVariant(final VcfHeader header, final VcfRecord rec) {
    super.tallyVariant(header, rec);

    // Handle the normal sample
    final int normalId = header.getSampleIndex(mNormalSampleName);
    if (isGenotype(rec, normalId, DIPLOID01)) {
      final Double rankingScore = rec.getSampleDouble(normalId, mRankingField);
      final String alleleDepth = rec.getSampleString(normalId, AD);
      insert(mNormalStore, rankingScore, alleleDepth);
    } else if (isGenotype(rec, normalId, DIPLOID00)) {
      // Handle the cancer sample
      final int cancerId = header.getSampleIndex(mCancerSampleName);
      if (isGenotype(rec, cancerId, DIPLOID01)) {
        final Double rankingScore = rec.getSampleDouble(cancerId, mRankingField);
        final String alleleDepth = rec.getSampleString(cancerId, AD);
        insert(mCancerStore, rankingScore, alleleDepth);
      }
    }
  }

  /*
   * In the ideal case, the number of alt alleles will be very close to the number of ref
   * alleles, in which care the correction factor will be 1.  However, in practice, the
   * number of alt alleles is less than the number of ref alleles and so the correction
   * is less than 1.  A correction greater than 1 which correspond to a strange case where
   * the alt allele was more common than the ref allele.
   */
  private double computeReferenceCorrection() {
    long ref = 0;
    long alt = 0;
    for (int k = 0; k < mNormalStore.size(); k++) {
      ref += mNormalStore.getRefCount(k);
      alt += mNormalStore.getAltCount(k);
    }
    final double total = ref + alt;
    final double correction = total == 0 ? 1 : 2.0 * alt / total;
    assert 0 <= correction;
    return correction;
  }

  private double computeContaminationEstimate() {
    long ref = 0;
    double alt = 0;
    for (int k = 0; k < mCancerStore.size(); k++) {
      ref += mCancerStore.getRefCount(k);
      alt += mCancerStore.getAltCount(k);
    }
    final double refBiasCorrection = computeReferenceCorrection();
    // If the reference bias is 0 (i.e. no alt evidence in the normal), then the best we can
    // do is to assume that there is no reference bias correction.
    if (refBiasCorrection > 0) {
      alt /= refBiasCorrection;
    }
    // A simple estimation based on 0/0 -> 0/1 calls in the derived sample, with a reference
    // bias correction applied to the alt allele count.  The final contamination is estimated
    // as 1 - 2 a / (a + r) where "a" is the corrected alt count and "r" is the ref count.
    final double total = ref + alt;
    // max below is a safety measure to handle the unlikely scenario where the alt count
    // exceeds the ref count.  Also, in the absence of any information we fallback to an
    // estimate of 0 contamination.
    final double contamEstimate = total == 0 ? 0 : Math.max(0, 1 - 2 * alt / total);
    assert 0 <= contamEstimate && contamEstimate <= 1 : contamEstimate + " " + alt + " " + ref;
    return contamEstimate;
  }

  @Override
  public String getStatistics() {
    final StringBuilder out = new StringBuilder();
    out.append("Passed Filters               : ").append(mTotalVariants).append(StringUtils.LS);
    out.append("Failed Filters               : ").append(mTotalFiltered).append(StringUtils.LS);
    if (mExcessCoverage > 0) {
      out.append("Excessive Coverage           : ").append(mExcessCoverage).append(StringUtils.LS);
    }
    if (mExcessHypotheses > 0) {
      out.append("Excessive Hypotheses         : ").append(mExcessHypotheses).append(StringUtils.LS);
    }
    if (mNoHypotheses > 0) {
      out.append("No Hypotheses                : ").append(mNoHypotheses).append(StringUtils.LS);
    }
    if (mCancerStore.size() != 0) {
      out.append("Estimated reference bias     : ").append(FORMAT_REAL.format(1 - computeReferenceCorrection())).append(StringUtils.LS);
      out.append("Estimated contamination      : ").append(FORMAT_REAL.format(computeContaminationEstimate())).append(StringUtils.LS);
    }

    final PerSampleVariantStatistics derivedSampleStats = mPerSampleStats.get(mCancerSampleName);
    if (derivedSampleStats != null) {
      derivedSampleStats.appendStatistics(out);
    }
    return out.toString();
  }
}
