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

package com.rtg.variant.bayes.multisample.cancer;

import java.util.Arrays;

import com.rtg.relation.Relationship;
import com.rtg.relation.Relationship.RelationshipType;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.VariantParams;
import com.rtg.variant.format.VcfFormatField;
import com.rtg.vcf.PerSampleVariantStatistics;
import com.rtg.vcf.VariantStatistics;
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

  private final String mRankingField;
  private final String mNormalSampleName;
  private final String mCancerSampleName;
  private final Scorer mNormalStore;
  private final Scorer mCancerStore;

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
    mNormalSampleName = genomeNames[AbstractSomaticCaller.NORMAL];
    mCancerSampleName = genomeNames[AbstractSomaticCaller.CANCER];
    mNormalStore = new TotalScorer();
    mCancerStore = new TopScorer(params.somaticParams().contaminationBasis());
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

  private void insert(final Scorer store, final Double rankingScore, final String alleleDepth) {
    final String[] alleleCounts = StringUtils.split(alleleDepth, ',');
    store.add(rankingScore, Integer.parseInt(alleleCounts[0]), Integer.parseInt(alleleCounts[1]));
  }

  /**
   * Increment statistics relevant to records even if they are not necessarily output.
   * @param header VCF header
   * @param rec the VCF record
   */
  public void countVariant(final VcfHeader header, final VcfRecord rec) {
    // Handle the normal sample
    final int normalId = header.getSampleIndex(mNormalSampleName);
    if (normalId != -1 && isGenotype(rec, normalId, DIPLOID01)) {
      final Double rankingScore = rec.getSampleDouble(normalId, mRankingField);
      final String alleleDepth = rec.getSampleString(normalId, AD);
      insert(mNormalStore, rankingScore, alleleDepth);
    } else if (normalId == -1 || isGenotype(rec, normalId, DIPLOID00)) {
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
   * In the ideal case, the count of alt allele evidence will be very close to the count of ref
   * allele evidence, in which case the alt allele fraction will be 0.5.  However, in practice, due
   * to reference bias during mapping the count for alt alleles is typically less than the count
   * for ref alleles and so the resulting alt fraction is less than 0.5.  A fraction greater than 0.5
   * corresponds to a strange case where the alt allele was more prevalent than the ref allele.
   */
  private double computeAltAlleleFraction() {
    final long ref = mNormalStore.getTotalRefCount();
    final long alt = mNormalStore.getTotalAltCount();
    final double total = ref + alt;
    final double altFraction = total == 0 ? 0.5 : alt / total;
    Diagnostic.developerLog("alt allele fraction: " + altFraction + " size=" + mNormalStore.size() + " altCount=" + alt + " refCount=" + ref);
    assert 0 <= altFraction && altFraction <= 1;
    return altFraction;
  }

  private double computeContaminationEstimate() {
    final long ref = mCancerStore.getTotalRefCount();
    final double alt = mCancerStore.getTotalAltCount();
    // A simple estimation based on 0/0 -> 0/1 calls in the derived sample, with a reference
    // bias correction applied to the alt allele count.  The final contamination is estimated
    // as 1 - 2 a / (a + r) where "a" is the corrected alt count and "r" is the ref count.
    final double total = ref + alt;
    // In the absence of any information assume the sample is pure.
    final double purityEstimate = total == 0 ? 1 : alt / total;
    final double normalAltAlleleFraction = computeAltAlleleFraction();
    // In the case where there is no alt allele information just take the raw purity
    // estimate.  Otherwise apply a correction based on the normal alt allele information.
    final double correctedPurityEstimate = normalAltAlleleFraction == 0 ? purityEstimate : Math.max(0, Math.min(1, purityEstimate / normalAltAlleleFraction));
    Diagnostic.developerLog("Purity estimate: " + correctedPurityEstimate + " size=" + mCancerStore.size() + " altCount=" + alt + " refCount=" + ref + " normalAltFraction=" + normalAltAlleleFraction);
    assert 0 <= correctedPurityEstimate && correctedPurityEstimate <= 1;
    return 1 - correctedPurityEstimate;
  }

  @Override
  public String getStatistics() {
    final StringBuilder out = new StringBuilder();
    out.append("Passed Filters               : ").append(mTotalPassed).append(StringUtils.LS);
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
      out.append("Alt fraction in 0/1 calls    : ").append(Utils.realFormat(computeAltAlleleFraction(), 2)).append(StringUtils.LS);
      out.append("Estimated contamination      : ").append(Utils.realFormat(computeContaminationEstimate(), 2)).append(StringUtils.LS);
    }

    final PerSampleVariantStatistics derivedSampleStats = mPerSampleStats.get(mCancerSampleName);
    if (derivedSampleStats != null) {
      derivedSampleStats.appendStatistics(out);
    }
    return out.toString();
  }
}
