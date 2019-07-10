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
package com.rtg.variant.format;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import com.rtg.reference.Ploidy;
import com.rtg.util.MathUtils;
import com.rtg.util.PosteriorUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.AlleleStatistics;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Statistics;
import com.rtg.variant.util.VariantUtils;
import com.rtg.vcf.VcfAnnotator;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.annotation.DerivedAnnotations;
import com.rtg.vcf.annotation.SplitContraryObservationAnnotator;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

/**
 * Enum of VCF FORMAT field implementations.
 */
public enum VcfFormatField {
  /**
   * Genotype.
   * Note: This Format field should always be populated first for each sample, and always be present as it also updates the ALT fields.
   * That means that this enum field should always be first.
   */
  GT(MetaType.STRING, VcfNumber.ONE, "Genotype") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      if (sample != null) {
        final String name = sample.getName();
        final Ploidy ploidy = sample.getPloidy();
        final Character previousNt = includePrevNt ? call.getLocus().getPreviousRefNt() : null;
        final int[] gts;
        if (name != null && !call.isFiltered(VariantFilter.FAILED_COMPLEX)) {
          final String ref = call.getLocus().getRefNts();
          if (sample.isIdentity()) {
            gts = ploidy == Ploidy.HAPLOID ? HAPLOID_REF_GT : DIPLOID_REF_GT;
          } else {
            final String[] cats = StringUtils.split(name, VariantUtils.COLON);
            gts = new int[cats.length];
            for (int k = 0; k < cats.length; ++k) {
              gts[k] = addAltAllele(cats[k], ref, previousNt, rec);
            }
            Arrays.sort(gts);
          }
        } else {
          gts = ploidy == Ploidy.HAPLOID ? HAPLOID_MISSING_GT : DIPLOID_MISSING_GT;
        }
        rec.addFormatAndSample(name(), VcfUtils.joinGt(false, gts));
      } else {
        rec.addFormatAndSample(name(), VcfRecord.MISSING);
      }
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return true;
    }
  },
  /**
   * Somatic variant allele. For cases where there is an allele which has a count or variant allelic frequency that exceeds call triggering threshold
   * but which isn't represented in the GT, we may want to add the allele in order to have it included in additional statistics such as <code>AD/VAC/VAF</code>.
   */
  VA(MetaType.INTEGER, VcfNumber.ONE, "Variant Allele") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      if (sample != null) {
        final String name = sample.getVariantAllele();
        if (name != null && !call.isFiltered(VariantFilter.FAILED_COMPLEX)) {
          final String ref = call.getLocus().getRefNts();
          final Character previousNt = includePrevNt ? call.getLocus().getPreviousRefNt() : null;
          final int gtnum = addAltAllele(name, ref, previousNt, rec);
          rec.addFormatAndSample(name(), String.valueOf(gtnum));
        }
      }
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getVariantAllele() != null;
    }
  },
  /** Coverage depth. */
  DP(MetaType.INTEGER, VcfNumber.ONE, "Read Depth") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      rec.addFormatAndSample(name(), String.valueOf(sample.getCoverage()));
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getCoverage() != null;
    }
  },
  /** Deviation from expected coverage. */
  DPR(MetaType.FLOAT, VcfNumber.ONE, "Ratio of read depth to expected read depth") {
    @Override
    protected void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      final double expected = params.expectedCoverage().expectedCoverage(call.getLocus().getSequenceName(), sampleName);
      rec.addFormatAndSample(name(), Utils.realFormat(sample.getCoverage() / expected, 3));
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getCoverage() != null && sampleName != null && params != null && params.expectedCoverage() != null && params.expectedCoverage().expectedCoverage(call.getLocus().getSequenceName(), sampleName) > 0;
    }
  },
  /** RTG error. */
  RE(MetaType.FLOAT, VcfNumber.ONE, "RTG Total Error") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      rec.addFormatAndSample(name(), Utils.realFormat(sample.getCorrection(), 3));
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getCorrection() != null;
    }
  },
  /** Ambiguity ratio. */
  AR(MetaType.FLOAT, VcfNumber.ONE, "Ambiguity Ratio") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      rec.addFormatAndSample(name(), Utils.realFormat(sample.getAmbiguityRatio(), 3));
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getAmbiguityRatio() != null;
    }
  },
  /** RTG sample quality. */
  RQ(MetaType.FLOAT, VcfNumber.ONE, "RTG sample quality") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      rec.addFormatAndSample(name(), Utils.realFormat(PosteriorUtils.phredIfy(sample.getNonIdentityPosterior()), 1));
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getNonIdentityPosterior() != null;
    }
  },
  /** Genotype quality. */
  GQ(MetaType.INTEGER, VcfNumber.ONE, "Genotype Quality") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      rec.addFormatAndSample(name(), MathUtils.cappedInt(MathUtils.round(PosteriorUtils.phredIfy(sample.getPosterior()))));
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getName() != null && !call.isFiltered(VariantFilter.FAILED_COMPLEX);
    }
  },
  /** RTG posterior. */
  RP(MetaType.FLOAT, VcfNumber.ONE, "RTG Posterior") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      rec.addFormatAndSample(name(), VariantUtils.formatPosterior(sample.getPosterior()));
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getName() != null && !call.isFiltered(VariantFilter.FAILED_COMPLEX);
    }
  },
  /** De novo allele. */
  DN(MetaType.CHARACTER, VcfNumber.ONE, "Indicates whether call is a putative de novo mutation") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      rec.addFormatAndSample(name(), sample.isDeNovo() == VariantSample.DeNovoStatus.IS_DE_NOVO ? "Y" : "N");
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return !call.isFiltered(VariantFilter.FAILED_COMPLEX) && sample != null && sample.getName() != null && sample.isDeNovo() != VariantSample.DeNovoStatus.UNSPECIFIED;
    }
  },
  /** De novo allele probability. */
  DNP(MetaType.FLOAT, VcfNumber.ONE, "Phred scaled probability that the call is due to a de novo mutation") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      rec.addFormatAndSample(name(), MathUtils.cappedInt(MathUtils.round(PosteriorUtils.phredIfy(sample.getDeNovoPosterior()))));
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.isDeNovo() != VariantSample.DeNovoStatus.UNSPECIFIED && sample.getDeNovoPosterior() != null;
    }
  },
  /** Hoeffding allele balance. */
  ABP(MetaType.FLOAT, VcfNumber.ONE, "Phred scaled probability that allele imbalance is present") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      rec.addFormatAndSample(name(), String.format(Locale.ROOT, "%3.2f", sample.getHoeffdingAlleleBalanceHom() != null ? sample.getHoeffdingAlleleBalanceHom() : sample.getHoeffdingAlleleBalanceHet()));
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && (sample.getHoeffdingAlleleBalanceHet() != null || sample.getHoeffdingAlleleBalanceHom() != null);
    }
  },
  /** Hoeffding strand bias. */
  SBP(MetaType.FLOAT, VcfNumber.ONE, "Phred scaled probability that strand bias is present") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      double score = sample.getHoeffdingStrandBiasAllele1();
      if (sample.getHoeffdingStrandBiasAllele2() != null && sample.getHoeffdingStrandBiasAllele2() > score) {
        score = sample.getHoeffdingStrandBiasAllele2();
      }
      rec.addFormatAndSample(name(), String.format(Locale.ROOT, "%3.2f", score));
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getHoeffdingStrandBiasAllele1() != null;
    }
  },
  /** Hoeffding read position bias. */
  RPB(MetaType.FLOAT, VcfNumber.ONE, "Phred scaled probability that read position bias is present") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      rec.addFormatAndSample(name(), String.format(Locale.ROOT, "%3.2f", sample.getHoeffdingReadPositionBias()));
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getHoeffdingReadPositionBias() != null;
    }
  },
  /** Hoeffding unmated bias. */
  PPB(MetaType.FLOAT, VcfNumber.ONE, "Phred scaled probability that there is a bias in the proportion of alignments that are properly paired") {
    @Override
    protected void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      double score = sample.getHoeffdingUnmatedBiasAllele1();
      if (sample.getHoeffdingUnmatedBiasAllele2() != null && sample.getHoeffdingUnmatedBiasAllele2() >  score) {
        score = sample.getHoeffdingUnmatedBiasAllele2();
      }
      rec.addFormatAndSample(name(), String.format(Locale.ROOT, "%3.2f", score));
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getHoeffdingUnmatedBiasAllele1() != null;
    }
  },
  /** Sum of quality for alleles evidence. */
  AQ(MetaType.FLOAT, VcfNumber.DOT, "Sum of quality for the evidence of the allele") {
    @Override
    protected void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      final String ref = includePrevNt ? rec.getRefCall().substring(1) : rec.getRefCall();
      final AlleleStatistics<?> counts = sample.getStats().counts();
      final Description desc = counts.getDescription();
      final int refHyp = desc.indexOf(ref);

      final StringBuilder value = new StringBuilder();
      value.append(refHyp == -1 ? 0 : Utils.realFormat(counts.qa(refHyp), 3));
      for (String altCall : rec.getAltCalls()) {
        final String adjAltCall = includePrevNt ? altCall.substring(1) : altCall;
        value.append(",").append(Utils.realFormat(counts.qa(desc.indexOf(adjAltCall)), 3));
      }
      rec.addFormatAndSample(name(), value.toString());
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getStats() != null;
    }
  },
  /** Placed unmapped ratio. */
  PUR(MetaType.FLOAT, VcfNumber.ONE, "Ratio of placed unmapped reads to mapped reads") {
    @Override
    protected void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      final double score = sample.getPlacedUnmappedRatio();
      rec.addFormatAndSample(name(), String.format(Locale.ROOT, "%3.2f", score));
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getPlacedUnmappedRatio() != null;
    }
  },
  /** RTG support statistics. */
  RS(MetaType.STRING, VcfNumber.DOT, "RTG Support Statistics") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      //Trim is to remove the leading TAB character
      rec.addFormatAndSample(name(), sample.getStatisticsString().trim().replaceAll("\t", ","));
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return !call.isComplexScored() && sample != null && sample.getStatisticsString() != null && sample.getStatisticsString().trim().length() > 0;
    }
  },
  /** Allelic depth, error-corrected. */
  ADE(MetaType.FLOAT, VcfNumber.DOT, "Allelic depths for the ref and alt alleles in the order listed, error corrected") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      final Statistics<?> stats = sample.getStats();
      final AlleleStatistics<?> counts = stats.counts();
      final Description description = counts.getDescription();
      final StringBuilder sb = new StringBuilder();
      final int refDescriptionCode = description.indexOf(getRefWithoutAnchor(rec, includePrevNt));
      sb.append(refDescriptionCode == -1 ? 0 : Utils.realFormat(counts.count(refDescriptionCode) - counts.error(refDescriptionCode), 1));
      for (final String altCall : rec.getAltCalls()) {
        final int altDescriptionCode = description.indexOf(includePrevNt ? altCall.substring(1) : altCall);
        sb.append(",").append(Utils.realFormat(counts.count(altDescriptionCode) - counts.error(altDescriptionCode), 1));
      }
      rec.addFormatAndSample(name(), sb.toString());
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getStats() != null;
    }
  },
  /** Allelic depth. */
  AD(MetaType.INTEGER, VcfNumber.DOT, "Allelic depths for the ref and alt alleles in the order listed") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      final Statistics<?> stats = sample.getStats();
      final Description description = stats.counts().getDescription();
      final StringBuilder sb = new StringBuilder();
      final int refDescriptionCode = description.indexOf(getRefWithoutAnchor(rec, includePrevNt));
      sb.append(refDescriptionCode == -1 ? 0 : MathUtils.round(stats.counts().count(refDescriptionCode)));
      for (final String altCall : rec.getAltCalls()) {
        final int altDescriptionCode = description.indexOf(includePrevNt ? altCall.substring(1) : altCall);
        sb.append(",").append(MathUtils.round(stats.counts().count(altDescriptionCode)));
      }
      rec.addFormatAndSample(name(), sb.toString());
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getStats() != null;
    }
  },
  /** Forward allelic depth. */
  ADF(MetaType.INTEGER, VcfNumber.DOT, "Allelic depth of reads on the forward strand") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      final Statistics<?> stats = sample.getStats();
      final AlleleStatistics<?> counts = stats.counts();
      final Description description = counts.getDescription();
      final StringBuilder sb = new StringBuilder();
      final int refDescriptionCode = description.indexOf(getRefWithoutAnchor(rec, includePrevNt));
      sb.append(refDescriptionCode == -1 ? 0 : Math.round(counts.forward(refDescriptionCode)));
      for (final String altCall : rec.getAltCalls()) {
        final int altDescriptionCode = description.indexOf(includePrevNt ? altCall.substring(1) : altCall);
        sb.append(",").append(Math.round(counts.forward(altDescriptionCode)));
      }
      rec.addFormatAndSample(name(), sb.toString());
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getStats() != null;
    }
  },
  /** Reverse allelic depth. */
  ADR(MetaType.INTEGER, VcfNumber.DOT, "Allelic depth of reads on the reverse strand") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      final Statistics<?> stats = sample.getStats();
      final AlleleStatistics<?> counts = stats.counts();
      final Description description = counts.getDescription();
      final StringBuilder sb = new StringBuilder();
      final int refDescriptionCode = description.indexOf(getRefWithoutAnchor(rec, includePrevNt));
      sb.append(refDescriptionCode == -1 ? 0 : Math.round(counts.backward(refDescriptionCode)));
      for (final String altCall : rec.getAltCalls()) {
        final int altDescriptionCode = description.indexOf(includePrevNt ? altCall.substring(1) : altCall);
        sb.append(",").append(Math.round(counts.backward(altDescriptionCode)));
      }
      rec.addFormatAndSample(name(), sb.toString());
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getStats() != null;
    }
  },
  /** Forward allelic depth for R1 reads. */
  ADF1(MetaType.INTEGER, VcfNumber.DOT, "Allelic depth of arm 1 reads on the forward strand") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      final Statistics<?> stats = sample.getStats();
      final AlleleStatistics<?> counts = stats.counts();
      final Description description = counts.getDescription();
      final StringBuilder sb = new StringBuilder();
      final int refDescriptionCode = description.indexOf(getRefWithoutAnchor(rec, includePrevNt));
      sb.append(refDescriptionCode == -1 ? 0 : Math.round(counts.forward1(refDescriptionCode)));
      for (final String altCall : rec.getAltCalls()) {
        final int altDescriptionCode = description.indexOf(includePrevNt ? altCall.substring(1) : altCall);
        sb.append(",").append(Math.round(counts.forward1(altDescriptionCode)));
      }
      rec.addFormatAndSample(name(), sb.toString());
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getStats() != null;
    }
  },
  /** Forward allelic depth for R2 reads. */
  ADF2(MetaType.INTEGER, VcfNumber.DOT, "Allelic depth of arm 2 reads on the forward strand") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      final Statistics<?> stats = sample.getStats();
      final AlleleStatistics<?> counts = stats.counts();
      final Description description = counts.getDescription();
      final StringBuilder sb = new StringBuilder();
      final int refDescriptionCode = description.indexOf(getRefWithoutAnchor(rec, includePrevNt));
      sb.append(refDescriptionCode == -1 ? 0 : Math.round(counts.forward2(refDescriptionCode)));
      for (final String altCall : rec.getAltCalls()) {
        final int altDescriptionCode = description.indexOf(includePrevNt ? altCall.substring(1) : altCall);
        sb.append(",").append(Math.round(counts.forward2(altDescriptionCode)));
      }
      rec.addFormatAndSample(name(), sb.toString());
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getStats() != null;
    }
  },
  /** Reverse allelic depth for arm 1 reads. */
  ADR1(MetaType.INTEGER, VcfNumber.DOT, "Allelic depth of arm 1 reads on the reverse strand") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      final Statistics<?> stats = sample.getStats();
      final AlleleStatistics<?> counts = stats.counts();
      final Description description = counts.getDescription();
      final StringBuilder sb = new StringBuilder();
      final int refDescriptionCode = description.indexOf(getRefWithoutAnchor(rec, includePrevNt));
      sb.append(refDescriptionCode == -1 ? 0 : Math.round(counts.backward1(refDescriptionCode)));
      for (final String altCall : rec.getAltCalls()) {
        final int altDescriptionCode = description.indexOf(includePrevNt ? altCall.substring(1) : altCall);
        sb.append(",").append(Math.round(counts.backward1(altDescriptionCode)));
      }
      rec.addFormatAndSample(name(), sb.toString());
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getStats() != null;
    }
  },
  /** Reverse allelic depth for R2 reads. */
  ADR2(MetaType.INTEGER, VcfNumber.DOT, "Allelic depth of arm 2 reads on the reverse strand") {
    @Override
    public void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      final Statistics<?> stats = sample.getStats();
      final AlleleStatistics<?> counts = stats.counts();
      final Description description = counts.getDescription();
      final StringBuilder sb = new StringBuilder();
      final int refDescriptionCode = description.indexOf(getRefWithoutAnchor(rec, includePrevNt));
      sb.append(refDescriptionCode == -1 ? 0 : Math.round(counts.backward2(refDescriptionCode)));
      for (final String altCall : rec.getAltCalls()) {
        final int altDescriptionCode = description.indexOf(includePrevNt ? altCall.substring(1) : altCall);
        sb.append(",").append(Math.round(counts.backward2(altDescriptionCode)));
      }
      rec.addFormatAndSample(name(), sb.toString());
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getStats() != null;
    }
  },
  /** Somatic score. */
  SSC(MetaType.FLOAT, VcfNumber.ONE, "Somatic score") {
    // Note TCGA VCF 1.2 spec defines this to be an integer 0-255, but that is pretty weak in terms of resolution
    @Override
    protected void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      if (sample.getDeNovoPosterior() != null) {
        rec.addFormatAndSample(name(), VariantUtils.formatPosterior(sample.getDeNovoPosterior()));
      }
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.isDeNovo() != VariantSample.DeNovoStatus.UNSPECIFIED && sample.getDeNovoPosterior() != null;
    }
  },
  /** Somatic status (loosely after TCGA "specification"). */
  SS(MetaType.INTEGER, VcfNumber.ONE, "Somatic status relative to original sample") {
    @Override
    protected void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      if (sample.isDeNovo() == VariantSample.DeNovoStatus.IS_DE_NOVO) {
        rec.addFormatAndSample(name(), "2");
      } else if (sample.isIdentity()) {
        rec.addFormatAndSample(name(), "0");
      } else {
        rec.addFormatAndSample(name(), "1");
      }
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.isDeNovo() != VariantSample.DeNovoStatus.UNSPECIFIED;
    }
  },

  // NOTE: VcfAnnotators (derived attributes) below here, non-derived above (convention is to have derived fields after non-derived)

  /** Variant allelic fraction. */
  VAF(DerivedAnnotations.VAF.getAnnotation()) {
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return rec.getFormat(ADE.name()) != null || rec.getFormat(AD.name()) != null;
    }
  },
  /** Genotype likelihood field (see VCF spec).  */
  GL(MetaType.FLOAT, VcfNumber.GENOTYPES, "Log_10 scaled genotype likelihoods. As defined in VCF specifications") {
    // GL is not a derived field but it depends upon the values of all sample's GT fields so must wait for the ALT arrays to be
    @Override
    protected void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
      final List<String> alleles = new ArrayList<>();
      alleles.add(rec.getRefCall());
      alleles.addAll(rec.getAltCalls());
      if (includePrevNt) {
        for (int i = 0; i < alleles.size(); ++i) {
          alleles.set(i, alleles.get(i).substring(1));
        }
      }
      final Map<Set<String>, Double> genotypeLikelihoods = sample.getGenotypeLikelihoods();
      if (genotypeLikelihoods != null && genotypeLikelihoods.size() > 0 && (sample.getPloidy() == Ploidy.HAPLOID || sample.getPloidy() == Ploidy.DIPLOID)) {
        final double[] likelihoods;
        likelihoods = sample.getPloidy() == Ploidy.DIPLOID ? GenotypeLikelihoodUtils.diploidLikelihoods(alleles, genotypeLikelihoods) : GenotypeLikelihoodUtils.haploidLikelihoods(alleles, genotypeLikelihoods);
        if (likelihoods == null) {
          Diagnostic.developerLog("Could not compute GL for variant at " + rec.getSequenceName() + ":" + (rec.getStart() + 1) + "-" + (rec.getEnd() + 1));
        } else {
          final StringBuilder sb = new StringBuilder();
          for (final double d : likelihoods) {
            if (sb.length() != 0) {
              sb.append(",");
            }
            sb.append(Utils.realFormat(d, 2));
          }
          rec.addFormatAndSample(name(), sb.toString());
        }
      }
    }
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null && sample.getGenotypeLikelihoods() != null && !sample.getGenotypeLikelihoods().isEmpty();
    }
    @Override
    public boolean isVcfAnnotator() {
      return true;
    }
  },
  /** GQ / DP. */
  GQD(DerivedAnnotations.GQD.getAnnotation()) {
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      final List<String> genotypeQuals = rec.getFormat(GQ.name());
      final List<String> depths = rec.getFormat(DP.name());
      if (genotypeQuals != null && depths != null && genotypeQuals.size() == depths.size()) {
        for (int i = 0; i < depths.size(); ++i) {
          if (!VcfRecord.MISSING.equals(genotypeQuals.get(i)) && !VcfRecord.MISSING.equals(depths.get(i))) {
            return true;
          }
        }
      }
      return false;
    }
  },
  /** Zygosity. */
  ZY(DerivedAnnotations.ZY.getAnnotation()) {
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      final List<String> genotypes = rec.getFormat(GT.name());
      if (genotypes != null) {
        for (final String gt : genotypes) {
          if (!gt.contains(VcfRecord.MISSING)) {
            return true;
          }
        }
      }
      return false;
    }
  },
  /** Ploidy. */
  PD(DerivedAnnotations.PD.getAnnotation()) {
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      final List<String> genotypes = rec.getFormat(GT.name());
      if (genotypes != null) {
        for (final String gt : genotypes) {
          if (!gt.contains(VcfRecord.MISSING)) {
            return true;
          }
        }
      }
      return false;
    }
  },
  /** Contrary observation count and fraction. */
  SCONT(new SplitContraryObservationAnnotator()) {
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      List<?> fld = rec.getFormat(AD.name());
      if (fld == null) {
        return false;
      }
      fld = rec.getFormat(GT.name());
      if (fld == null) {
        return false;
      }
      fld = rec.getFormat(SS.name());
      if (fld != null) {
        return true;
      }
      fld = rec.getFormat(DN.name());
      if (fld != null) {
        return true;
      }
      return false;
    }
  },
  /** Sum of quality of the alternate observations. */
  QA(DerivedAnnotations.QA.getAnnotation()) {
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null; // todo this condition looks bogus
    }
  },
  /** Difference in mean quality for called alleles. */
  MEANQAD(DerivedAnnotations.MEANQAD.getAnnotation()) {
    @Override
    public boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params) {
      return sample != null; // todo this condition looks bogus
    }
  }
  ;

  private static String getRefWithoutAnchor(final VcfRecord rec, final boolean includePrevNt) {
    final String ref = rec.getRefCall();
    return includePrevNt ? ref.substring(1) : ref;
  }

  private static final int[] HAPLOID_REF_GT = {0};
  private static final int[] DIPLOID_REF_GT = {0, 0};
  private static final int[] HAPLOID_MISSING_GT = {VcfUtils.MISSING_GT};
  private static final int[] DIPLOID_MISSING_GT = {VcfUtils.MISSING_GT, VcfUtils.MISSING_GT};

  private final VcfAnnotator mAnnotator;
  private final MetaType mMetaType;
  private final VcfNumber mVcfNumber;
  private final String mDescription;

  VcfFormatField(final VcfAnnotator annotator) {
    mAnnotator = annotator;
    mMetaType = null;
    mVcfNumber = null;
    mDescription = null;
  }

  VcfFormatField(final MetaType type, final VcfNumber number, final String description) {
    mAnnotator = null;
    mMetaType = type;
    mVcfNumber = number;
    mDescription = description;
  }

  /**
   * Update the VCF header with the field description.
   * @param header the VCF header for which the field description will be added.
   */
  public void updateHeader(final VcfHeader header) {
    if (mAnnotator != null) {
      mAnnotator.updateHeader(header);
    } else {
      header.addFormatField(name(), mMetaType, mVcfNumber, mDescription);
    }
  }

  /**
   * Update the VCF record with the value for the FORMAT field.
   * Should only be called if at least one sample for the call has a value for the given FORMAT field.
   * @param rec the VCF record for which the FORMAT field value will be added.
   * @param call the variant data to use to update the VCF record.
   * @param sample the variant sample to output the value for.
   * @param sampleName the name of the sample to update.
   * @param params the variant output options params.
   * @param includePrevNt true when including previous NT in VCF output, false otherwise.
   */
  protected void updateRecordSample(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
    //One location for most of the MISSING record output.
    if (hasValue(rec, call, sample, sampleName, params)) {
      updateVcfRecord(rec, call, sample, sampleName, params, includePrevNt);
    } else {
      rec.addFormatAndSample(name(), VcfRecord.MISSING);
    }
  }

  /**
   * Update the VCF record with the values for the FORMAT field.
   * Should only be called if at least one sample for the call the given FORMAT field.
   * @param rec the VCF record for which the FORMAT field value will be added.
   * @param call the variant data to use to update the VCF record.
   * @param sampleNames the list of sample names in the output VCF.
   * @param params the variant output options params.
   * @param includePrevNt true when including previous NT in VCF output, false otherwise.
   */
  public void updateRecord(VcfRecord rec, Variant call, String[] sampleNames, VariantParams params, boolean includePrevNt) {
    if (isVcfAnnotator() && hasValue(rec, call, null, null, params)) {
      updateRecordSample(rec, call, null, null, params, includePrevNt);
    } else {
      //Should process all samples for each field one field at a time
      for (int i = 0; i < sampleNames.length; ++i) {
        updateRecordSample(rec, call, call.getSample(i), sampleNames[i], params, includePrevNt);
      }
    }
  }

  /**
   * Update the VCF record with the value for the FORMAT field.
   * It can be assumed that the record has output to print.
   * @param rec the VCF record for which the FORMAT field value will be added.
   * @param call the variant data to use to update the VCF record.
   * @param sample the variant sample to output the value for.
   * @param sampleName the name of the sample to update.
   * @param params the variant output options params.
   * @param includePrevNt true when including previous NT in VCF output, false otherwise.
   */
  protected void updateVcfRecord(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params, boolean includePrevNt) {
    mAnnotator.annotate(rec);
  }

  /**
   * Find if this sample has information for format field.
   * @param rec the VCF record being annotated
   * @param call the call the sample belongs to.
   * @param sample the sample to check.
   * @param sampleName the name of the sample to check.
   * @param params the variant output options params.
   * @return true if there is information in the sample for the format field.
   */
  public abstract boolean hasValue(VcfRecord rec, Variant call, VariantSample sample, String sampleName, VariantParams params);

  /**
   * Check if this enum value uses a VCF annotator and should be run after all normal
   * parts of the VCF have been constructed.
   * @return true if this enum value uses information from the VCF record to annotate, false if it uses the Variant instead
   */
  public boolean isVcfAnnotator() {
    return mAnnotator != null;
  }

  protected static int addAltAllele(String alt, String ref, Character previousNt, VcfRecord rec) {
    if (alt.equals(ref)) {
      return 0;
    } else {
      final StringBuilder alternate = new StringBuilder();
      if (previousNt != null) {
        alternate.append(previousNt);
      }
      alternate.append(alt);
      int pos = rec.getAltCalls().indexOf(alternate.toString());
      if (pos < 0) {
        pos = rec.addAltCall(alternate.toString()).getAltCalls().size() - 1;
      }
      return pos + 1;
    }
  }
}
