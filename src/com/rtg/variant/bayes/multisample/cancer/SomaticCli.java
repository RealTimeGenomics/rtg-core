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

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ParamsTask;
import com.rtg.reference.Sex;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.Relationship.RelationshipType;
import com.rtg.usage.UsageMetric;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.avr.AbstractPredictModel;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;
import com.rtg.variant.bayes.multisample.MultisampleTask;
import com.rtg.variant.format.VcfFormatField;

/**
 * Somatic caller
 */
public class SomaticCli extends AbstractMultisampleCli {

  /** Module name displayed on help */
  public static final String MODULE_NAME = "somatic";

  private static final String SOMATIC_FLAG = "somatic";
  private static final String SOMATIC_PRIORS_FLAG = "somatic-priors";
  private static final String LOH_FLAG = "loh";
  // Either this ...
  private static final String PEDIGREE_FLAG = "Xpedigree";
  // or these
  private static final String DERIVED_FLAG = "derived";
  private static final String ORIGINAL_FLAG = "original";
  private static final String CONTAMINATION_FLAG = "contamination";

  private static final String SEX_FLAG = "sex";
  private static final String INCLUDE_GERMLINE_FLAG = "Xinclude-germline";
  private static final String INCLUDE_GAIN_OF_REFERENCE = "include-gain-of-reference";

  private static final String MISMATCHED_PARAMS_ERROR1 = "Cannot use --" + PEDIGREE_FLAG + " in conjunction with --" + DERIVED_FLAG + ", --" + ORIGINAL_FLAG + ", or --" + CONTAMINATION_FLAG;
  private static final String MISMATCHED_PARAMS_ERROR2 = "All of --" + DERIVED_FLAG + ", --" + ORIGINAL_FLAG + ", and --" + CONTAMINATION_FLAG + " must be specified";

  static class MultigenomeValidator implements Validator {
    @Override
    public boolean isValid(final CFlags flags) {
      if (!AbstractMultisampleCli.validateCommonOptions(flags)) {
        return false;
      }
      if (flags.isSet(PEDIGREE_FLAG)) {
        if (flags.isSet(DERIVED_FLAG) || flags.isSet(ORIGINAL_FLAG) || flags.isSet(CONTAMINATION_FLAG)) {
          flags.error(MISMATCHED_PARAMS_ERROR1);
          return false;
        }
      } else {
        if (!flags.isSet(DERIVED_FLAG) || !flags.isSet(ORIGINAL_FLAG) || !flags.isSet(CONTAMINATION_FLAG)) {
          flags.error(MISMATCHED_PARAMS_ERROR2);
          return false;
        }
        final String original = (String) flags.getValue(ORIGINAL_FLAG);
        final String derived = (String) flags.getValue(DERIVED_FLAG);
        if (original.equals(derived)) {
          flags.error("Original and derived must be different samples");
          return false;
        }
      }
      if (Sex.EITHER != flags.getValue(SEX_FLAG) && !CommonFlags.validateSexTemplateReference(flags, SEX_FLAG, null, TEMPLATE_FLAG)) {
        return false;
      }
      final double som = (Double) flags.getValue(SOMATIC_FLAG);
      if (som <= 0 || som >= 1) {
        flags.error("--" + SOMATIC_FLAG + " should be a probability 0<s<1");
        return false;
      }
      final double loh = (Double) flags.getValue(LOH_FLAG);
      if (loh < 0 || loh > 1) {
        flags.error("--" + LOH_FLAG + " should be a probability 0<=s<=1");
        return false;
      }
      return true;
    }
  }

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  protected void initFlags() {
    initLocalFlags(mFlags);
  }

  @Override
  protected double defaultCoverageCutoffMultiplier() {
    return 25;
  }

  void initLocalFlags(CFlags flags) {
    initFlags(flags);
    CommonFlags.initAvrModel(flags, false);
    CommonFlags.initMinAvrScore(flags);
    flags.setDescription("Performs a somatic variant analysis.");
    flags.setValidator(new MultigenomeValidator());
    final Flag inFlag = flags.registerRequired(File.class, "file", "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    final Flag listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, "FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);
    flags.registerOptional('r', PEDIGREE_FLAG, File.class, "file", "relationships file").setCategory(INPUT_OUTPUT);
    final Flag contamFlag = flags.registerOptional(CONTAMINATION_FLAG, Double.class, "float", "estimated fraction of contamination in derived sample").setCategory(SENSITIVITY_TUNING);
    final Flag derivedFlag = flags.registerOptional(DERIVED_FLAG, String.class, "string", "sample identifier used in read groups for derived sample").setCategory(INPUT_OUTPUT);
    final Flag originalFlag = flags.registerOptional(ORIGINAL_FLAG, String.class, "string", "sample identifier used in read groups for original sample").setCategory(INPUT_OUTPUT);
    flags.registerOptional(SEX_FLAG, Sex.class, "sex", "sex of individual", Sex.EITHER).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('s', SOMATIC_FLAG, Double.class, "float", "default prior probability of a somatic SNP mutation in the derived sample", 1e-6).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(SOMATIC_PRIORS_FLAG, File.class, "file", "if set, use the BED file to generated site specific somatic priors").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(LOH_FLAG, Double.class, "float", "prior probability that a loss of heterozygosity event has occurred", 0.0).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(INCLUDE_GERMLINE_FLAG, "include germline variants in output VCF").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('G', INCLUDE_GAIN_OF_REFERENCE, "include gain of reference somatic calls in output VCF").setCategory(SENSITIVITY_TUNING);
    AbstractMultisampleCli.registerComplexPruningFlags(flags);
    flags.addRequiredSet(derivedFlag, originalFlag, contamFlag, inFlag);
    flags.addRequiredSet(derivedFlag, originalFlag, contamFlag, listFlag);
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
      ped.addRelationship(RelationshipType.ORIGINAL_DERIVED, original, derived).setProperty(CONTAMINATION_FLAG, String.valueOf(mFlags.getValue(CONTAMINATION_FLAG)));
      return ped;
    }
  }

  @Override
  protected VariantParamsBuilder makeParamsBuilder() throws InvalidParamsException, IOException {
    final VariantParamsBuilder vpb = super.makeParamsBuilder();
    vpb
      .somaticRate((Double) mFlags.getValue(SOMATIC_FLAG))
      .interestingThreshold(0)
      .includeGermlineVariants(mFlags.isSet(INCLUDE_GERMLINE_FLAG))
      .includeGainOfReference(mFlags.isSet(INCLUDE_GAIN_OF_REFERENCE))
      .lohPrior((Double) mFlags.getValue(LOH_FLAG))
      .sex((Sex) mFlags.getValue(SEX_FLAG));
    if (mFlags.isSet(SOMATIC_PRIORS_FLAG)) {
      final PriorBedRangeLoader loader = new PriorBedRangeLoader();
      final File sspBedFile = (File) mFlags.getValue(SOMATIC_PRIORS_FLAG);
      loader.loadRanges(getSimpleRegionRestriction(), sspBedFile);
      final ReferenceRanges<Double> ssp = loader.getReferenceRanges();
      vpb.siteSpecificSomaticPriors(ssp);
      Diagnostic.userLog("Loaded site specific somatic priors from " + sspBedFile);
    }
    return vpb;
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
