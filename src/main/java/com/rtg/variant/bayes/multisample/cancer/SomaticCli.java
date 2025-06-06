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

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ParamsTask;
import com.rtg.reference.Sex;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.Relationship;
import com.rtg.relation.Relationship.RelationshipType;
import com.rtg.usage.UsageMetric;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.SomaticParamsBuilder;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.avr.AbstractPredictModel;
import com.rtg.variant.avr.AvrUtils;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;
import com.rtg.variant.bayes.multisample.MultisampleTask;
import com.rtg.variant.format.VcfFormatField;

/**
 * Somatic caller entry point.
 */
public class SomaticCli extends AbstractMultisampleCli {

  // The name of the default AVR model used for somatic calling
  static final String SOMATIC_MODEL_DEFAULT = "illumina-somatic.avr";

  protected static final String SOMATIC_FLAG = "somatic";
  private static final String SOMATIC_PRIORS_FLAG = "somatic-priors";
  protected static final String LOH_FLAG = "loh";
  // Either this ...
  protected static final String PEDIGREE_FLAG = "Xpedigree";
  // or these
  protected static final String DERIVED_FLAG = "derived";
  // Unregistered for tumor only caller so made package local
  static final String ORIGINAL_FLAG = "original";
  protected static final String CONTAMINATION_FLAG = Relationship.CONTAMINATION;

  protected static final String REVERSE_CONTAMINATION_FLAG = "X" + Relationship.REVERSE_CONTAMINATION;

  protected static final String SEX_FLAG = "sex";
  static final String INCLUDE_GERMLINE_FLAG = "include-germline";
  private static final String INCLUDE_GAIN_OF_REFERENCE = "include-gain-of-reference";
  private static final String USE_SOMATIC_ALLELIC_FRACTION = "enable-somatic-allelic-fraction";
  private static final String CONTAMINATION_BASIS_FLAG = "Xcontamination-basis";

  private static final String MISMATCHED_PARAMS_ERROR1 = "Cannot use --" + PEDIGREE_FLAG + " in conjunction with --" + DERIVED_FLAG + ", --" + ORIGINAL_FLAG + ", or --" + CONTAMINATION_FLAG;
  private static final String MISMATCHED_PARAMS_ERROR2 = "All of --" + DERIVED_FLAG + ", --" + ORIGINAL_FLAG + ", and --" + CONTAMINATION_FLAG + " must be specified";
  private static final int DEFAULT_CONTAMINATION_BASIS = 1000;

  private static class MultigenomeValidator implements Validator {
    @Override
    public boolean isValid(final CFlags flags) {
      if (!AbstractMultisampleCli.validateCommonOptions(flags)) {
        return false;
      }
      if (flags.isSet(PEDIGREE_FLAG)) {
        if (flags.isSet(DERIVED_FLAG) || flags.isSet(ORIGINAL_FLAG) || flags.isSet(CONTAMINATION_FLAG)) {
          flags.setParseMessage(MISMATCHED_PARAMS_ERROR1);
          return false;
        }
      } else {
        if (!flags.isSet(DERIVED_FLAG) || !flags.isSet(ORIGINAL_FLAG) || !flags.isSet(CONTAMINATION_FLAG)) {
          flags.setParseMessage(MISMATCHED_PARAMS_ERROR2);
          return false;
        }
        final String original = (String) flags.getValue(ORIGINAL_FLAG);
        final String derived = (String) flags.getValue(DERIVED_FLAG);
        if (original.equals(derived)) {
          flags.setParseMessage("Original and derived must be different samples");
          return false;
        }
      }
      if (Sex.EITHER != flags.getValue(SEX_FLAG) && !CommonFlags.validateSexTemplateReference(flags, SEX_FLAG, null, CommonFlags.TEMPLATE_FLAG)) {
        return false;
      }
      if (!flags.checkInRange(SOMATIC_FLAG, 0.0, false, 1.0, false)
        || !flags.checkInRange(LOH_FLAG, 0.0, true, 1.0, true)
        || !flags.checkInRange(CONTAMINATION_FLAG, 0.0, true, 1.0, false)
        || !flags.checkInRange(REVERSE_CONTAMINATION_FLAG, 0.0, true, 1.0, false)) {
        return false;
      }
      return true;
    }
  }

  @Override
  public String moduleName() {
    return "somatic";
  }

  @Override
  public String description() {
    return "call variants for a tumor/normal pair";
  }

  @Override
  protected void initFlags() {
    initLocalFlags(mFlags);
  }

  @Override
  protected double defaultCoverageCutoffMultiplier() {
    return 25;
  }
  void commonSomaticFlags(CFlags flags) {
    flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, CommonFlags.FILE, "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);
    flags.registerOptional(CONTAMINATION_FLAG, Double.class, CommonFlags.FLOAT, "estimated fraction of contamination in derived sample").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(REVERSE_CONTAMINATION_FLAG, Double.class, CommonFlags.FLOAT, "estimated fraction of derived sample in original sample", 0.0).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(SEX_FLAG, Sex.class, "sex", "sex of individual", Sex.EITHER).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('s', SOMATIC_FLAG, Double.class, CommonFlags.FLOAT, "default prior probability of a somatic SNP mutation in the derived sample", 1e-6).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(SOMATIC_PRIORS_FLAG, File.class, CommonFlags.FILE, "if set, use the BED file to generate site specific somatic priors").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(LOH_FLAG, Double.class, CommonFlags.FLOAT, "prior probability that a loss of heterozygosity event has occurred", 0.0).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('G', INCLUDE_GAIN_OF_REFERENCE, "include gain of reference somatic calls in output VCF").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(USE_SOMATIC_ALLELIC_FRACTION, "if set, incorporate the expected somatic allelic fraction in scoring").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(CONTAMINATION_BASIS_FLAG, Integer.class, CommonFlags.INT, "number of examples to use in computing the contamination estimate", DEFAULT_CONTAMINATION_BASIS).setCategory(SENSITIVITY_TUNING);
    final Flag<File> inFlag = flags.registerRequired(File.class, CommonFlags.FILE, "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    registerAllelicTriggers(flags);
    registerComplexPruningFlags(flags, true);
  }

  void initLocalFlags(CFlags flags) {
    initFlags(flags);
    AvrUtils.initAvrModel(flags, false, SOMATIC_MODEL_DEFAULT);
    CommonFlags.initMinAvrScore(flags);
    commonSomaticFlags(flags);
    flags.registerOptional(DERIVED_FLAG, String.class, CommonFlags.STRING, "sample identifier used in read groups for derived sample").setCategory(INPUT_OUTPUT);
    flags.registerOptional(ORIGINAL_FLAG, String.class, CommonFlags.STRING, "sample identifier used in read groups for original sample").setCategory(INPUT_OUTPUT);
    flags.registerOptional(INCLUDE_GERMLINE_FLAG, "include germline variants in output VCF").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('r', PEDIGREE_FLAG, File.class, CommonFlags.FILE, "relationships file").setCategory(INPUT_OUTPUT);
    flags.setDescription("Performs a somatic variant analysis.");
    flags.setValidator(new MultigenomeValidator());
    requiredSet(flags);
  }

  private void requiredSet(CFlags flags) {
    flags.addRequiredSet(flags.getFlag(DERIVED_FLAG), flags.getFlag(ORIGINAL_FLAG), flags.getFlag(CONTAMINATION_FLAG), flags.getAnonymousFlag(0));
    flags.addRequiredSet(flags.getFlag(DERIVED_FLAG), flags.getFlag(ORIGINAL_FLAG), flags.getFlag(CONTAMINATION_FLAG), flags.getFlag(CommonFlags.INPUT_LIST_FLAG));
  }

  @Override
  protected GenomeRelationships grf() throws IOException {
    if (mFlags.isSet(PEDIGREE_FLAG)) {
      return GenomeRelationships.loadGenomeRelationships((File) mFlags.getValue(PEDIGREE_FLAG));
    } else {
      assert mFlags.isSet(ORIGINAL_FLAG);
      assert mFlags.isSet(DERIVED_FLAG);
      assert mFlags.isSet(CONTAMINATION_FLAG);

      final GenomeRelationships ped = new GenomeRelationships();
      final String original = (String) mFlags.getValue(ORIGINAL_FLAG);
      final String derived = (String) mFlags.getValue(DERIVED_FLAG);
      final Sex sex = (Sex) mFlags.getValue(SEX_FLAG);
      ped.addGenome(original, sex);
      ped.addGenome(derived, sex);
      final Relationship relationship = ped.addRelationship(RelationshipType.ORIGINAL_DERIVED, original, derived);
      relationship.setProperty(Relationship.CONTAMINATION, String.valueOf(mFlags.getValue(CONTAMINATION_FLAG)));
      relationship.setProperty(Relationship.REVERSE_CONTAMINATION, String.valueOf(mFlags.getValue(REVERSE_CONTAMINATION_FLAG)));
      return ped;
    }
  }

  @Override
  protected VariantParamsBuilder makeParamsBuilder() throws IOException {
    final SomaticParamsBuilder somaticBuilder = getSomaticParamsBuilder();
    final VariantParamsBuilder vpb = super.makeParamsBuilder();
    vpb
      .interestingThreshold(0)
      .genomePriors(GenomePriorParams.builder().contraryProbability((Double) mFlags.getValue(X_CONTRARY_FLAG)).create())
      .sex((Sex) mFlags.getValue(SEX_FLAG));
    if (mFlags.isSet(SOMATIC_PRIORS_FLAG)) {
      final PriorBedRangeLoader loader = new PriorBedRangeLoader();
      final File sspBedFile = (File) mFlags.getValue(SOMATIC_PRIORS_FLAG);
      loader.loadRanges(getSimpleRegionRestriction(), sspBedFile);
      final ReferenceRanges<Double> ssp = loader.getReferenceRanges();
      somaticBuilder.siteSpecificSomaticPriors(ssp);
      Diagnostic.userLog("Loaded site specific somatic priors from " + sspBedFile);
    }
    vpb.somaticParams(somaticBuilder.create());
    return vpb;
  }

  protected SomaticParamsBuilder getSomaticParamsBuilder() {
    final SomaticParamsBuilder somaticBuilder = new SomaticParamsBuilder();
    somaticBuilder
      .somaticRate((Double) mFlags.getValue(SOMATIC_FLAG))
      .includeGermlineVariants(mFlags.isSet(INCLUDE_GERMLINE_FLAG))
      .includeGainOfReference(mFlags.isSet(INCLUDE_GAIN_OF_REFERENCE))
      .somaticAlleleBalance(mFlags.isSet(USE_SOMATIC_ALLELIC_FRACTION))
      .contaminationBasis((Integer) mFlags.getValue(CONTAMINATION_BASIS_FLAG))
      .lohPrior((Double) mFlags.getValue(LOH_FLAG));
    return somaticBuilder;
  }

  @Override
  protected SomaticStatistics getStatistics(VariantParams params) {
    return new SomaticStatistics(params, params.avrModelFile() != null ? AbstractPredictModel.AVR : VcfFormatField.GQ.name());
  }

  @Override
  public ParamsTask<?, ?> task(final VariantParams params, final OutputStream out) throws IOException {
    final UsageMetric usageMetric = mUsageMetric == null ? new UsageMetric() : mUsageMetric; //create when null to cover some testing
    return new MultisampleTask<>(params, new SomaticCallerConfiguration.Configurator(), out, getStatistics(params), usageMetric);
  }
}
