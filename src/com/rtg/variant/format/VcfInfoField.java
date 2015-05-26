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

import com.rtg.util.Pair;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.variant.CoverageThreshold;
import com.rtg.variant.PosteriorUtils;
import com.rtg.variant.StaticThreshold;
import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.VariantParams;
import com.rtg.variant.util.VariantUtils;
import com.rtg.vcf.VcfAnnotator;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.annotation.DerivedAnnotations;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

/**
 * Enum of VCF INFO field implementations
 */
public enum VcfInfoField {
  /** Somatic Mutation */
  SOMATIC {
    @Override
    public void updateHeader(VcfHeader header) {
      header.addInfoField(name(), MetaType.STRING, new VcfNumber("1"), "Indicates the variant is a somatic mutation");
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      if (call.getNormalCancerScore() != null) {
        rec.addInfo(name(), formatPossibleCause(call, includePrevNt));
      }
    }
  },
  /** Loss Of Heterozygosity */
  LOH {
    @Override
    public void updateHeader(VcfHeader header) {
      header.addInfoField(name(), MetaType.FLOAT, new VcfNumber("1"), "Indicates whether or not variant is a potential loss of heterozygosity");
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      if (call.getLoh() != null) {
        rec.addInfo(name(), Utils.realFormat(call.getLoh()));
      }
    }
  },
  /** RTG Normal Cancer Score */
  NCS {
    @Override
    public void updateHeader(VcfHeader header) {
      header.addInfoField(name(), MetaType.FLOAT, new VcfNumber("1"), "RTG normal cancer score");
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      if (call.getNormalCancerScore() != null) {
        rec.addInfo(name(), Utils.realFormat(PosteriorUtils.phredIfy(Math.abs(call.getNormalCancerScore())), 1));
      }
    }
  },
  /** Disease Mutation */
  DISEASE {
    @Override
    public void updateHeader(VcfHeader header) {
      header.addInfoField(name(), MetaType.STRING, new VcfNumber("1"), "Indicates the variant is linked to the disease");
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      if (call.getDiseasePresenceScore() != null) {
        rec.addInfo(name(), formatPossibleCause(call, includePrevNt));
      }
    }
  },
  /** RTG Disease Score */
  RDS {
    @Override
    public void updateHeader(VcfHeader header) {
      header.addInfoField(name(), MetaType.FLOAT, new VcfNumber("1"), "RTG disease call score");
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      if (call.getPossibleCauseScore() != null) {
        rec.addInfo(name(), VariantUtils.formatPosterior(call.getPossibleCauseScore()));
      }
    }
  },
  /** RTG Disease Presence Score */
  DPS {
    @Override
    public void updateHeader(VcfHeader header) {
      header.addInfoField(name(), MetaType.FLOAT, new VcfNumber("1"), "RTG disease presence score");
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      if (call.getDiseasePresenceScore() != null) {
        rec.addInfo(name(), Utils.realFormat(PosteriorUtils.phredIfy(Math.abs(call.getDiseasePresenceScore())), 1));
      }
    }
  },
  /** Combined Depth */
  DP {
    @Override
    public void updateHeader(VcfHeader header) {
      header.addInfoField(name(), MetaType.INTEGER, new VcfNumber("1"), "Combined read depth for variant over all samples");
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      final Pair<Integer, Boolean> cov = getCoverage(call);
      if (cov.getB()) {
        rec.addInfo(name(), Integer.toString(cov.getA()));
      }
    }
  },
  /** Deviation from expected combined coverage */
  DPR {
    @Override
    public void updateHeader(VcfHeader header) {
      header.addInfoField(name(), MetaType.FLOAT, new VcfNumber("1"), "Ratio of combined read depth for variant to expected combined read depth");
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      final Pair<Integer, Boolean> cov = getCoverage(call);
      if (cov.getB() && params != null && params.expectedCoverage() != null) {
        final double expected = params.expectedCoverage().expectedTotalCoverage(call.getLocus().getSequenceName());
        if (expected > 0) {
          rec.addInfo(name(), Utils.realFormat(cov.getA() / expected, 3));
        }
      }
    }
  },
  /** Complex Called */
  XRX {
    @Override
    public void updateHeader(VcfHeader header) {
      header.addInfoField(name(), MetaType.FLAG, new VcfNumber("0"), "RTG variant was called using complex caller");
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      if (call.isComplexScored()) {
        rec.addInfo(name());
      }
    }
  },
  /** Equivalent Call */
  RCE {
    @Override
    public void updateHeader(VcfHeader header) {
      header.addInfoField(name(), MetaType.FLAG, new VcfNumber("0"), "RTG variant is equivalent to the previous variant");
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      if (call.isComplexEquivalent()) {
        rec.addInfo(name());
      }
    }
  },
  /** Coverage Threshold */
  CT {
    @Override
    public void updateHeader(VcfHeader header) {
      header.addInfoField(name(), MetaType.INTEGER, new VcfNumber("1"), "Coverage threshold that was applied");
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      if (call.isFiltered(VariantFilter.COVERAGE)) {
        final CoverageThreshold threshold = params != null ? params.maxCoverageFilter() : DUMMY_COVERAGE_THRESHOLD;
        rec.addInfo(name(), Integer.toString(threshold.thresholdSingle(call.getLocus().getSequenceName())));
      }
    }
  },
  /** Call Trimmed */
  RTRM {
    @Override
    public void updateHeader(VcfHeader header) {
      header.addInfoField(name(), MetaType.FLAG, new VcfNumber("0"), "Complex called variant was trimmed");
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      if (call.isTrimmed()) {
        rec.addInfo(name());
      }
    }
  },
  /** Call Split */
  RSPLT {
    @Override
    public void updateHeader(VcfHeader header) {
      header.addInfoField(name(), MetaType.INTEGER, new VcfNumber("1"), "Complex called variant was split");
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      if (call.getSplitId() > 0) {
        rec.addInfo(name(), Integer.toString(call.getSplitId()));
      }
    }
  },
  /** No reference hypothesis */
  NREF {
    @Override
    public void updateHeader(VcfHeader header) {
      header.addInfoField(name(), MetaType.FLAG, new VcfNumber("0"), "REF was not considered a valid allele during calling");
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      if (call.isInvalidRef()) {
        rec.addInfo(name());
      }
    }
  },

  // NOTE: VcfAnnotators (derived attributes) below here, non-derived above (convention is to have derived fields after non-derived)

  /** Inbreeding Coefficient */
  IC {
    @Override
    public void updateHeader(VcfHeader header) {
      IC_ANNOTATOR.updateHeader(header);
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      IC_ANNOTATOR.annotate(rec);
    }
    @Override
    public boolean isVcfAnnotator() {
      return true;
    }
  },
  /** Hardy-Weinberg Equilibrium Probability */
  EP {
    @Override
    public void updateHeader(VcfHeader header) {
      EP_ANNOTATOR.updateHeader(header);
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      EP_ANNOTATOR.annotate(rec);
    }
    @Override
    public boolean isVcfAnnotator() {
      return true;
    }
  },
  /** Length of the longest allele */
  LAL {
    @Override
    public void updateHeader(VcfHeader header) {
      LAL_ANNOTATOR.updateHeader(header);
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      LAL_ANNOTATOR.annotate(rec);
    }
    @Override
    public boolean isVcfAnnotator() {
      return true;
    }
  },
  /** QUAL / DP */
  QD {
    @Override
    public void updateHeader(VcfHeader header) {
      QD_ANNOTATOR.updateHeader(header);
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      QD_ANNOTATOR.annotate(rec);
    }
    @Override
    public boolean isVcfAnnotator() {
      return true;
    }
  },
  /** Number of alternative alleles */
  NAA {
    @Override
    public void updateHeader(VcfHeader header) {
      NAA_ANNOTATOR.updateHeader(header);
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      NAA_ANNOTATOR.annotate(rec);
    }
    @Override
    public boolean isVcfAnnotator() {
      return true;
    }
  },
  /** Allele count in genotypes, for each alternative allele */
  AC {
    @Override
    public void updateHeader(VcfHeader header) {
      AC_ANNOTATOR.updateHeader(header);
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      AC_ANNOTATOR.annotate(rec);
    }
    @Override
    public boolean isVcfAnnotator() {
      return true;
    }
  },
  /** Total number of alleles in called genotypes */
  AN {
    @Override
    public void updateHeader(VcfHeader header) {
      AN_ANNOTATOR.updateHeader(header);
    }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) {
      AN_ANNOTATOR.annotate(rec);
    }
    @Override
    public boolean isVcfAnnotator() {
      return true;
    }
  },
  /** Segregation probability annotator (a flag to turn on, no implementation) */
  SGP {
    //Do nothing, will be turned on and used by family caller using this as a flag setting
    @Override
    public void updateHeader(VcfHeader header) { }
    @Override
    public void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt) { }
    @Override
    public boolean isVcfAnnotator() {
      return true;
    }
  };

  private static final CoverageThreshold DUMMY_COVERAGE_THRESHOLD = new StaticThreshold(0);
  private static final VcfAnnotator IC_ANNOTATOR = VcfUtils.getAnnotator(DerivedAnnotations.IC);
  private static final VcfAnnotator EP_ANNOTATOR = VcfUtils.getAnnotator(DerivedAnnotations.EP);
  private static final VcfAnnotator LAL_ANNOTATOR = VcfUtils.getAnnotator(DerivedAnnotations.LAL);
  private static final VcfAnnotator QD_ANNOTATOR = VcfUtils.getAnnotator(DerivedAnnotations.QD);
  private static final VcfAnnotator NAA_ANNOTATOR = VcfUtils.getAnnotator(DerivedAnnotations.NAA);
  private static final VcfAnnotator AC_ANNOTATOR = VcfUtils.getAnnotator(DerivedAnnotations.AC);
  private static final VcfAnnotator AN_ANNOTATOR = VcfUtils.getAnnotator(DerivedAnnotations.AN);

  /**
   * Update the VCF header with the field description.
   * @param header the VCF header for which the field description will be added.
   */
  public abstract void updateHeader(VcfHeader header);

  /**
   * Update the VCF record with the value for the INFO field.
   * @param rec the VCF record for which the INFO field value will be added.
   * @param call the variant data to use to update the VCF record.
   * @param params the variant output options params.
   * @param includePrevNt true when including previous NT in VCF output, false otherwise.
   */
  public abstract void updateRecord(VcfRecord rec, Variant call, VariantParams params, boolean includePrevNt);

  /**
   * Check if this enum value uses a VCF annotator and should be run after all normal
   * parts of the VCF have been constructed.
   * @return true if this enum value uses information from the VCF record to annotate, false if it uses the Variant instead
   */
  public boolean isVcfAnnotator() {
    return false;
  }

  protected static String formatPossibleCause(Variant call, boolean includePreviousNt) {
    if (call.getPossibleCause() == null) {
      return "*";
    }
    if (includePreviousNt) {
      final String[] parts = StringUtils.split(call.getPossibleCause(), ':');
      if (parts.length == 2) {
        final char nt = call.getLocus().getPreviousRefNt();
        return nt + parts[0] + ":" + nt + parts[1];
      }
      return call.getLocus().getPreviousRefNt() + call.getPossibleCause();
    }
    return call.getPossibleCause();
  }

  /**
   * Get the coverage information from the variant.
   * @param call the variant call to get the coverage information from.
   * @return a pair object with the A part holding the coverage and the B part holding if there were coverage values.
   */
  protected static Pair<Integer, Boolean> getCoverage(Variant call) {
    int coverage = 0;
    boolean hasCoverage = false;
    for (int i = 0; i < call.getNumberOfSamples(); i++) {
      if (call.getSample(i) != null && call.getSample(i).getCoverage() != null) {
        coverage += call.getSample(i).getCoverage();
        hasCoverage = true;
      }
    }
    return new Pair<>(coverage, hasCoverage);
  }
}
