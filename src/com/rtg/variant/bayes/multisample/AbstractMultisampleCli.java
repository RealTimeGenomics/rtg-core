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
package com.rtg.variant.bayes.multisample;


import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collection;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.calibrate.CalibratedPerSequenceExpectedCoverage;
import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.SamCalibrationInputs;
import com.rtg.launcher.BuildCommon;
import com.rtg.launcher.CommandLineFiles;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.OutputParams;
import com.rtg.launcher.ParamsCli;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.SdfUtils;
import com.rtg.reference.SexMemo;
import com.rtg.relation.GenomeRelationships;
import com.rtg.sam.SamFilterOptions;
import com.rtg.sam.SamFilterParams;
import com.rtg.sam.SamRangeUtils;
import com.rtg.sam.SamUtils;
import com.rtg.util.IORunnable;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.variant.CalibratedPerSequenceThreshold;
import com.rtg.variant.CalibratedPerSequenceThreshold.ThresholdFunction;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.StaticThreshold;
import com.rtg.variant.ThreadingEnvironment;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.avr.AvrUtils;
import com.rtg.variant.bayes.AlleleBalanceProbability;
import com.rtg.variant.bayes.BinomialAlleleBalance;
import com.rtg.variant.bayes.HoeffdingAlleleBalance;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.format.VcfFormatField;
import com.rtg.variant.format.VcfInfoField;
import com.rtg.vcf.VariantStatistics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;

/**
 * Common stuff for multi-genome SNP caller.
 */
@TestClass({"com.rtg.variant.bayes.multisample.family.FamilyCliTest", "com.rtg.variant.bayes.multisample.singleton.SingletonCliTest"})
public abstract class AbstractMultisampleCli extends ParamsCli<VariantParams> {

  /** Default absolute coverage cutoff. */
  private static final int DEFAULT_COVERAGE_CUTOFF = 200;
  private static final String ALL_FLAG = "all";
  private static final String MACHINE_ERRORS_FLAG = "machine-errors";

  // Flags for applying filters to output records.
  private static final String FILTER_AMBIGUITY_FLAG = "filter-ambiguity";
  private static final String FILTER_DEPTH_FLAG = "filter-depth";
  private static final String FILTER_DEPTH_MULTIPLIER_FLAG = "filter-depth-multiplier";
  private static final String BED_FILTER_FLAG = "filter-bed";
  static final String COVERAGE_BYPASS_FLAG = "max-coverage";
  private static final String COVERAGE_BYPASS_MULTIPLIER_FLAG = "max-coverage-multiplier";
  /** flag for disabling use of calibration files */
  public static final String NO_CALIBRATION = "no-calibration";
  private static final String R_DEFAULT_FLAG = "rdefault-mated";
  private static final String UNMATED_R_DEFAULT_FLAG = "rdefault-unmated";
  private static final String SNPS_ONLY_FLAG = "snps-only";
  private static final String POPULATION_PRIORS = "population-priors";
  protected static final String MIN_VARIANT_ALLELIC_DEPTH = "min-variant-allelic-depth";
  protected static final String MIN_VARIANT_ALLELIC_FRACTION = "min-variant-allelic-fraction";
  /** Minimum read base quality filter flag */
  public static final String MIN_BASE_QUALITY_FLAG = "min-base-quality";

  private static final String X_PRIORS_FLAG = "Xpriors";
  private static final String X_Q_DEFAULT_FLAG = "Xqdefault";
  private static final String X_INDEL_TRIGGER_FRACTION_FLAG = "Xindel-trigger-fraction";
  private static final String X_IGNORE_QUALITIES_FLAG = "Xignore-qualities";
  private static final String X_CHUNKING_FLAG = "Xchunking";
  private static final String X_LOOKAHEAD_FLAG = "Xlookahead";
  private static final String X_VCF_RP = "Xvcfrp";
  private static final String X_HYPER_COMPLEX_LENGTH_FLAG = "Xhyper-complex-length";
  private static final String X_INTERESTING_SEPARATION_FLAG = "Xinteresting-separation";
  private static final String X_SIMPLE_REPEAT_EXTENSION = "Xsimple-repeat-extension";
  private static final String X_R_IGNORE_FLAG = "Xrignore";
  private static final String X_R_MAX_FLAG = "Xrmax";
  private static final String X_UNMATED_R_MAX_FLAG = "Xunmated-rmax";
  private static final String X_THREADING_ENVIRONMENT = "Xthreading-env";
  private static final String X_IO_THREADS = "Xio-threads";
  protected static final String X_PRUNE_HYPOTHESES = "Xprune-hypotheses";
  protected static final String X_MAX_COMPLEX_HYPOTHESES = "Xmax-complex-hypotheses";
  private static final String X_NO_COMPLEX_CALLS_FLAG = "Xno-complex-calls";
  private static final String X_TRIM_SPLIT_FLAG = "Xtrim-split";
  private static final String X_ALT_MULTIPLIER_FLAG = "Xalt-coverage-multiplier";
  protected static final String X_CONTRARY_FLAG = "Xcontrary-probability";

  private static final String X_IGNORE_SAM_HEADER_INCOMPATIBILITY = "Xignore-incompatible-sam-headers";
  private static final String X_ALTERNATE_SAM_HEADER = "Xalternate-sam-header";

  private static final String X_INFO_ANNOTATION_FLAG = "Xinfo-annotation";
  private static final String X_FORMAT_ANNOTATION_FLAG = "Xformat-annotation";
  protected static final String USE_ALLELIC_FRACTION = "enable-allelic-fraction";
  private static final String X_ALLELE_BALANCE_PROBABILITY = "Xallele-balance-probability"; // Fine grained control over the method used
  private static final String X_EXPECTED_ALLELIC_FRACTION = "Xexpected-allele-balance"; // Fine grained control over the expected fraction
  private static final String X_EXPAND_COMPLEX_READ_QUERIES = "Xexpand-complex-read-queries";
  private static final String X_COMPLEX_USE_SOFT_CLIP = "Xcomplex-use-soft-clip";
  private static final String X_FORCE_COMPLEX_REGION = "Xforce-complex-region";
  private static final double LOW_COVERAGE_WARNING_THRESHOLD = 2.0;

  /**
   * validate common flags
   * @param flags flags to validate
   * @return true of validated false otherwise
   */
  public static boolean validateCommonOptions(final CFlags flags) {
    if (!CommonFlags.validateOutputDirectory(flags)) {
      return false;
    }
    if (!CommonFlags.validateTemplate(flags)) {
      return false;
    }
    if (!CommonFlags.checkFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, Integer.MAX_VALUE, true)) {
      return false;
    }
    if (!CommonFlags.validateThreads(flags)) {
      return false;
    }
    if (!SamFilterOptions.validateFilterFlags(flags, false)) {
      return false;
    }
    if (flags.isSet(FILTER_AMBIGUITY_FLAG)) {
      final IntegerOrPercentage amb = (IntegerOrPercentage) flags.getValue(FILTER_AMBIGUITY_FLAG);
      if (amb.getRawValue() < 0) {
        Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + FILTER_AMBIGUITY_FLAG, amb + "", "0%");
        return false;
      } else if (amb.getRawValue() > 100) {
        Diagnostic.error(ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, "--" + FILTER_AMBIGUITY_FLAG, amb + "", "100%");
        return false;
      }
    }
    if (flags.isSet(MIN_BASE_QUALITY_FLAG)) {
      final Integer minBaseQuality = (Integer) flags.getValue(MIN_BASE_QUALITY_FLAG);
      if (minBaseQuality < 0) {
        Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + MIN_BASE_QUALITY_FLAG,  minBaseQuality + "", "0");
        return false;
      }
    }
    // Chunks sizes less than twice the maximum read length can create problems in the
    // circular buffer used in the reading of SAM record.
    final Integer chunkSize = (Integer) flags.getValue(X_CHUNKING_FLAG);
    if (chunkSize < 1000) {
      flags.setParseMessage("Invalid chunk size, minimum chunking size is 1000");
      return false;
    }
    final Integer lookahead = (Integer) flags.getValue(X_LOOKAHEAD_FLAG);
    if (lookahead < 2) {
      flags.setParseMessage("Invalid lookahead value, minimum is 2");
      return false;
    }
    if (flags.isSet(X_THREADING_ENVIRONMENT)) {
      final String val = (String) flags.getValue(X_THREADING_ENVIRONMENT);
      if (!"simple".equalsIgnoreCase(val) && !"single".equalsIgnoreCase(val) && !"parallel".equalsIgnoreCase(val) && !val.toLowerCase(Locale.getDefault()).startsWith("random")) {
        flags.setParseMessage("Unknown value " + val + " for --" + X_THREADING_ENVIRONMENT);
        return false;
      }
      if (val.toLowerCase(Locale.getDefault()).startsWith("random")) {
        final String[] split = val.split("=");
        if (split.length != 2) {
          flags.setParseMessage("A valid seed is required for random threading environment");
          return false;
        }
        boolean valid = true;
        try {
          Long.parseLong(split[1]);
        } catch (final NumberFormatException ex) {
          valid = false;
        }
        if (!valid) {
          flags.setParseMessage("A valid seed is required for random threading environment");
          return false;
        }
      }
    }
    if (flags.isSet(FILTER_DEPTH_FLAG) && flags.isSet(FILTER_DEPTH_MULTIPLIER_FLAG)) {
      flags.setParseMessage("Only one of --" + FILTER_DEPTH_FLAG + " or --" + FILTER_DEPTH_MULTIPLIER_FLAG + " can be set");
      return false;
    }
    if (flags.isSet(COVERAGE_BYPASS_FLAG) && flags.isSet(COVERAGE_BYPASS_MULTIPLIER_FLAG)) {
      flags.setParseMessage("Only one of --" + COVERAGE_BYPASS_FLAG + " or --" + COVERAGE_BYPASS_MULTIPLIER_FLAG + " can be set");
      return false;
    }
    if ((flags.isSet(FILTER_DEPTH_FLAG) && flags.isSet(COVERAGE_BYPASS_MULTIPLIER_FLAG))
        || (flags.isSet(FILTER_DEPTH_MULTIPLIER_FLAG) && flags.isSet(COVERAGE_BYPASS_FLAG))) {
      flags.setParseMessage("Cannot mix ratio based and fixed coverage thresholding");
      return false;
    }
    if (flags.isSet(FILTER_DEPTH_FLAG)) {
      if ((Integer) flags.getValue(FILTER_DEPTH_FLAG) > (Integer) flags.getValue(COVERAGE_BYPASS_FLAG)) {
        Diagnostic.warning("--" + FILTER_DEPTH_FLAG + " is set higher than --" + COVERAGE_BYPASS_FLAG + ", consider increasing --" + COVERAGE_BYPASS_FLAG);
      }
    }
    if (flags.isSet(FILTER_DEPTH_MULTIPLIER_FLAG)) {
      if ((Double) flags.getValue(FILTER_DEPTH_MULTIPLIER_FLAG) > (Double) flags.getValue(COVERAGE_BYPASS_MULTIPLIER_FLAG)) {
        Diagnostic.warning("--" + FILTER_DEPTH_MULTIPLIER_FLAG + " is set higher than --" + COVERAGE_BYPASS_MULTIPLIER_FLAG + ", consider increasing --" + COVERAGE_BYPASS_MULTIPLIER_FLAG);
      }
    }
    if (flags.isSet(POPULATION_PRIORS)) {
      final File f = (File) flags.getValue(POPULATION_PRIORS);
      if (!f.exists()) {
        flags.setParseMessage("Given file for flag --" + POPULATION_PRIORS + " " + f.getPath() + " does not exist");
        return false;
      }
    }
    if (!flags.checkInRange(CommonFlags.FILTER_AVR_FLAG, 0.0, true, 1.0, true)
      || !flags.checkInRange(X_CONTRARY_FLAG, 0.0, false, 1.0, true)) {
      return false;
    }
    if (flags.isSet(BED_FILTER_FLAG)) {
      final File bedRegionsFile = (File) flags.getFlag(BED_FILTER_FLAG).getValue();
      if (!bedRegionsFile.exists()) {
        Diagnostic.error(ErrorType.FILE_NOT_FOUND,
            "The specified file, \"" + bedRegionsFile.getPath() + "\", does not exist.");
        return false;
      }
    }
    return true;
  }

  /**
   * Sam file flags (input/file list)
   * @param flags object to add flag to
   */
  public static void addSamFileFlags(CFlags flags) {
    final Flag<File> inFlag = flags.registerRequired(File.class, CommonFlags.FILE, "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    final Flag<File> listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, CommonFlags.FILE, "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);
    flags.addRequiredSet(inFlag);
    flags.addRequiredSet(listFlag);
  }

  protected double defaultCoverageCutoffMultiplier() {
    return 5;
  }

  /**
   * initialise flags. Don't add things here unless they are common to ALL of family, somatic, and disease callers.
   * @param flags flags to initialise
   */
  public void initFlags(CFlags flags) {
    CommonFlagCategories.setCategories(flags);
    flags.registerOptional(NO_CALIBRATION, "if set, ignore mapping calibration files").setCategory(UTILITY);
    CommonFlags.initReferenceTemplate(flags, true);
    CommonFlags.initOutputDirFlag(flags);
    CommonFlags.initNoGzip(flags);
    CommonFlags.initThreadsFlag(flags);
    CommonFlags.initIndexFlags(flags);
    SamFilterOptions.registerMinMapQFlag(flags);
    SamFilterOptions.registerRestrictionFlag(flags);
    SamFilterOptions.registerBedRestrictionFlag(flags);
    SamFilterOptions.registerKeepDuplicatesFlag(flags);
    flags.registerOptional('q', X_Q_DEFAULT_FLAG, Integer.class, CommonFlags.INT, "for reads that have no quality information use this as the default quality (in Phred format from 0 to 63)", 20).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('m', MACHINE_ERRORS_FLAG, String.class, CommonFlags.STRING, "if set, force sequencer machine settings. Allowed values are [default, illumina, ls454_se, ls454_pe, complete, iontorrent]").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(R_DEFAULT_FLAG, Integer.class, CommonFlags.INT, "for mated reads that have no mapping quality supplied use this as the default quality (in Phred format from 0 to 63)", 20).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(UNMATED_R_DEFAULT_FLAG, Integer.class, CommonFlags.INT, "for unmated reads that have no mapping quality supplied use this as the default quality (in Phred format from 0 to 63)", 20).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(POPULATION_PRIORS, File.class, CommonFlags.FILE, "if set, use the VCF file to generate population based site-specific priors").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(COVERAGE_BYPASS_FLAG, Integer.class, CommonFlags.INT, "skip calling in sites with per sample read depth exceeding this value", DEFAULT_COVERAGE_CUTOFF).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(COVERAGE_BYPASS_MULTIPLIER_FLAG, Double.class, CommonFlags.FLOAT, "skip calling in sites with combined depth exceeding multiplier * average combined coverage determined from calibration", defaultCoverageCutoffMultiplier()).setCategory(SENSITIVITY_TUNING); //set a coverage threshold for use with the average coverage determined by calibration files. The threshold will be the average coverage + (multiplier * square root of average coverage)
    flags.registerOptional(USE_ALLELIC_FRACTION, "if set, incorporate the expected allelic fraction in scoring").setCategory(SENSITIVITY_TUNING);

    flags.registerOptional(SNPS_ONLY_FLAG, "if set, will output simple SNPs only").setCategory(REPORTING);
    flags.registerOptional('a', ALL_FLAG, "write variant calls covering every position irrespective of thresholds").setCategory(REPORTING);
    flags.registerOptional(FILTER_AMBIGUITY_FLAG, IntegerOrPercentage.class, CommonFlags.INT, "threshold for ambiguity filter applied to output variants").setCategory(REPORTING);
    flags.registerOptional(FILTER_DEPTH_FLAG, Integer.class, CommonFlags.INT, "apply a fixed depth of coverage filter to output variants").setCategory(REPORTING);
    flags.registerOptional(FILTER_DEPTH_MULTIPLIER_FLAG, Double.class, CommonFlags.FLOAT, "apply a ratio based depth filter. The filter will be multiplier * average coverage determined from calibration files").setCategory(REPORTING); //set a coverage threshold for use with the average coverage determined by calibration files. The threshold will be the average coverage + (multiplier * square root of average coverage)
    flags.registerOptional(BED_FILTER_FLAG, File.class, CommonFlags.FILE, "apply a position based filter, retaining only variants that fall in these BED regions").setCategory(REPORTING);
    flags.registerOptional(MIN_BASE_QUALITY_FLAG, Integer.class, CommonFlags.INT, "phred scaled quality score, read bases below this quality will be treated as unknowns", 0).setCategory(SENSITIVITY_TUNING);

    flags.registerOptional(X_ALT_MULTIPLIER_FLAG, "determine coverage thresholds using avg_cov + multiplier * sqrt(avg_cov)").setCategory(REPORTING); /* old multiplier * avg_cov instead of as */

    flags.registerOptional(X_PRIORS_FLAG, String.class, CommonFlags.STRING, "selects a properties file specifying the priors. Either a file name or one of [human]", "human").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(X_IGNORE_QUALITIES_FLAG, "if set, will ignore quality scores associated with reads and use the default instead").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(X_INDEL_TRIGGER_FRACTION_FLAG, Double.class, CommonFlags.INT, "if set, fraction of evidence at a position that must be indel to trigger complex calling", 0.05).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(X_CHUNKING_FLAG, Integer.class, CommonFlags.INT, "number of nucleotide positions considered per chunk", 1000).setCategory(UTILITY);
    flags.registerOptional(X_LOOKAHEAD_FLAG, Integer.class, CommonFlags.INT, "number of chunks to prefetch", 2).setCategory(UTILITY);
    flags.registerOptional(X_VCF_RP, "include RTG posterior in VCF output").setCategory(REPORTING);
    flags.registerOptional(X_HYPER_COMPLEX_LENGTH_FLAG, Integer.class, CommonFlags.INT, "the length beyond which complex regions are considered hyper complex", Integer.MAX_VALUE).setCategory(SENSITIVITY_TUNING);
    // flags.registerOptional(X_INTERESTING_THRESHOLD_FLAG, Double.class, CommonFlags.FLOAT, "posterior threshold below which an identity call is considered interesting", Double.valueOf(1.0 / VariantUtils.LOG_10)).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(X_INTERESTING_SEPARATION_FLAG, Integer.class, CommonFlags.INT, "the maximum distance over which two interesting calls will be considered part of the same complex region", 4).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(X_SIMPLE_REPEAT_EXTENSION, Boolean.class, CommonFlags.BOOL, "extend complex regions over DNA simple repeats", Boolean.TRUE).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(X_R_IGNORE_FLAG, "if set any supplied read qualities will be ignored and the defaults used").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(X_R_MAX_FLAG, Integer.class, CommonFlags.INT, "for mated reads this is the maximum value that the read quality can have", 255).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(X_UNMATED_R_MAX_FLAG, Integer.class, CommonFlags.INT, "for unmated reads this is the maximum value that the read quality can have", 255).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(X_THREADING_ENVIRONMENT, String.class, CommonFlags.STRING, "threading environment to be used. One of [single, random=seed, parallel]", "parallel").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(X_IO_THREADS, Integer.class, CommonFlags.INT, "number of threads to use for IO (Default is the number of available cores)").setCategory(UTILITY);
    flags.registerOptional(X_NO_COMPLEX_CALLS_FLAG, "turn off attempting calls in complex region").setCategory(INPUT_OUTPUT);
    flags.registerOptional(X_TRIM_SPLIT_FLAG, DecomposerType.class, CommonFlags.STRING, "type of normalization/decomposition applied to haplotype calls", DecomposerType.ALIGN).setCategory(REPORTING);
    flags.registerOptional(X_IGNORE_SAM_HEADER_INCOMPATIBILITY, "ignore incompatible SAM headers when merging SAM results").setCategory(UTILITY);
    flags.registerOptional(X_ALTERNATE_SAM_HEADER, File.class, CommonFlags.FILE, "treat all SAM records as having the supplied header").setCategory(UTILITY);
    flags.registerOptional(X_CONTRARY_FLAG, Double.class, CommonFlags.FLOAT, "probability used to penalize contrary evidence in somatic calls", 0.01).setCategory(SENSITIVITY_TUNING);

    //Extra INFO / FORMAT fields
    flags.registerOptional(X_INFO_ANNOTATION_FLAG, VcfInfoField.class, CommonFlags.STRING, "additional VCF INFO fields").setCategory(REPORTING)
      .setParameterRange(new String[]{VcfInfoField.IC.name(), VcfInfoField.EP.name(), VcfInfoField.LAL.name(), VcfInfoField.QD.name(), VcfInfoField.NAA.name(), VcfInfoField.AN.name(), VcfInfoField.AC.name(), VcfInfoField.RTRM.name(), VcfInfoField.RSPLT.name(), VcfInfoField.SGP.name()}).setMaxCount(Integer.MAX_VALUE).enableCsv();
    flags.registerOptional(X_FORMAT_ANNOTATION_FLAG, VcfFormatField.class, CommonFlags.STRING, "additional VCF FORMAT fields").setCategory(REPORTING)
      .setParameterRange(new String[]{VcfFormatField.GQD.name(), VcfFormatField.ZY.name(), VcfFormatField.PD.name(), VcfFormatField.ADF.name(), VcfFormatField.ADR.name(), VcfFormatField.ADF1.name(), VcfFormatField.ADF2.name(), VcfFormatField.ADR1.name(), VcfFormatField.ADR2.name()}).setMaxCount(Integer.MAX_VALUE).enableCsv();

    flags.registerOptional(X_ALLELE_BALANCE_PROBABILITY, AlleleBalanceFactor.class, CommonFlags.STRING, "method for calculating allele balance", AlleleBalanceFactor.NONE).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(X_EXPECTED_ALLELIC_FRACTION, Double.class, CommonFlags.FLOAT, "expected germline heterozygous alternate allele fraction", 0.5).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(X_EXPAND_COMPLEX_READ_QUERIES, Boolean.class, CommonFlags.BOOL, "expand queries for reads by one base either side of a complex region", Boolean.FALSE).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(X_COMPLEX_USE_SOFT_CLIP, Boolean.class, CommonFlags.BOOL, "use soft clipped bases in evidence for complex calls", Boolean.TRUE).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(X_FORCE_COMPLEX_REGION, String.class, CommonFlags.STRING, "Force a complex region over specified range").setCategory(UTILITY);
  }

  protected void registerAllelicTriggers(CFlags flags) {
    flags.registerOptional(MIN_VARIANT_ALLELIC_DEPTH, Double.class, CommonFlags.FLOAT, "if set, also output sites that meet this minimum quality-adjusted alternate allelic depth").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(MIN_VARIANT_ALLELIC_FRACTION, Double.class, CommonFlags.FLOAT, "if set, also output sites that meet this minimum quality-adjusted alternate allelic fraction").setCategory(SENSITIVITY_TUNING);
  }

  /**
   * Add flags for complex hypothesis pruning to the specified flags object
   * @param flags flags object to modify
   * @param pruneDefault whether pruning should be enabled by default
   */
  protected void registerComplexPruningFlags(CFlags flags, boolean pruneDefault) {
    flags.registerOptional(X_PRUNE_HYPOTHESES, Boolean.class, CommonFlags.BOOL, "prune hypotheses during complex calling", pruneDefault).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(X_MAX_COMPLEX_HYPOTHESES, Integer.class, CommonFlags.INT, "maximum number of remaining hypotheses after pruning", 21).setCategory(SENSITIVITY_TUNING);
  }

  protected abstract GenomeRelationships grf() throws IOException;

  @Override
  public VariantParams makeParams() throws IOException {
    final VariantParams localParams = makeParamsBuilder().create();
    localParams.globalIntegrity();
    return localParams;
  }

  protected VariantParamsBuilder makeParamsBuilder() throws IOException {
    final VariantParamsBuilder builder = new VariantParamsBuilder();
    builder.name(mFlags.getName());
    builder.callLevel(mFlags.isSet(ALL_FLAG) ? VariantOutputLevel.ALL : VariantOutputLevel.INTERESTING);
    builder.chunkSize((Integer) mFlags.getValue(X_CHUNKING_FLAG));
    builder.lookAhead((Integer) mFlags.getValue(X_LOOKAHEAD_FLAG));
    builder.hyperComplexLength((Integer) mFlags.getValue(X_HYPER_COMPLEX_LENGTH_FLAG));
    builder.interestingSeparation((Integer) mFlags.getValue(X_INTERESTING_SEPARATION_FLAG));
    builder.simpleRepeatExtension((Boolean) mFlags.getValue(X_SIMPLE_REPEAT_EXTENSION));
    builder.ignoreReadQuality(mFlags.isSet(X_R_IGNORE_FLAG));
    builder.outputIndex(!mFlags.isSet(CommonFlags.NO_INDEX));
    builder.noComplexCalls(mFlags.isSet(X_NO_COMPLEX_CALLS_FLAG));
    builder.outputNonSnps(!mFlags.isSet(SNPS_ONLY_FLAG));
    builder.trimSplit((DecomposerType) mFlags.getValue(X_TRIM_SPLIT_FLAG));
    builder.vcfRp(mFlags.isSet(X_VCF_RP));
    builder.avrModelFile(AvrUtils.getAvrModel(mFlags, false));
    builder.ignoreQualityScores(mFlags.isSet(X_IGNORE_QUALITIES_FLAG));
    final SamFilterParams filterParams = SamFilterOptions.makeFilterParamsBuilder(mFlags).excludeUnplaced(true).minLength(1).excludeVariantInvalid(true).create();
    builder.filterParams(filterParams);
    if (mFlags.isSet(BED_FILTER_FLAG)) {
      builder.regionsFilterBedFile((File) mFlags.getValue(BED_FILTER_FLAG));
    }
    parseThreads(builder);
    processDefaultAndMaxMappingQuality(builder);
    final boolean zip = !mFlags.isSet(CommonFlags.NO_GZIP);
    final OutputParams outParams = new OutputParams((File) mFlags.getValue(CommonFlags.OUTPUT_FLAG), mFlags.isSet(BuildCommon.PROGRESS_FLAG), zip);
    builder.outputParams(outParams);
    if (mFlags.getFlag(X_PRUNE_HYPOTHESES) != null) {
      builder.pruneHypotheses((Boolean) mFlags.getValue(X_PRUNE_HYPOTHESES));
    }
    if (mFlags.getFlag(X_MAX_COMPLEX_HYPOTHESES) != null && mFlags.isSet(X_MAX_COMPLEX_HYPOTHESES)) {
      builder.maxComplexHypotheses((Integer) mFlags.getValue(X_MAX_COMPLEX_HYPOTHESES));
    }
    builder.indelTriggerFraction((Double) mFlags.getValue(X_INDEL_TRIGGER_FRACTION_FLAG));
    if (mFlags.isSet(FILTER_AMBIGUITY_FLAG)) {
      final IntegerOrPercentage amb = (IntegerOrPercentage) mFlags.getValue(FILTER_AMBIGUITY_FLAG);
      builder.maxAmbiguity(amb.getRawValue() / 100.0);
    }
    final EnumSet<VcfInfoField> infos = EnumSet.noneOf(VcfInfoField.class);
    for (Object obj : mFlags.getValues(X_INFO_ANNOTATION_FLAG)) {
      infos.add((VcfInfoField) obj);
    }
    final EnumSet<VcfFormatField> formats = EnumSet.noneOf(VcfFormatField.class);
    for (Object obj : mFlags.getValues(X_FORMAT_ANNOTATION_FLAG)) {
      formats.add((VcfFormatField) obj);
    }
    builder.infoAnnotations(infos).formatAnnotations(formats);
    if (mFlags.getFlag(CommonFlags.FILTER_AVR_FLAG) != null && mFlags.isSet(CommonFlags.FILTER_AVR_FLAG)) {
      builder.minAvrScore((Double) mFlags.getValue(CommonFlags.FILTER_AVR_FLAG));
    }
    if (mFlags.getFlag(MIN_BASE_QUALITY_FLAG) != null && mFlags.isSet(MIN_BASE_QUALITY_FLAG)) {
      builder.minBaseQuality((Integer) mFlags.getValue(MIN_BASE_QUALITY_FLAG));
    }
    builder.genomePriors(GenomePriorParams.builder().genomePriors((String) mFlags.getValue(X_PRIORS_FLAG)).contraryProbability((Double) mFlags.getValue(X_CONTRARY_FLAG)).create());
    if (mFlags.isSet(MACHINE_ERRORS_FLAG)) {
      builder.machineErrorName((String) mFlags.getValue(MACHINE_ERRORS_FLAG));
    }
    if (mFlags.isSet(POPULATION_PRIORS)) {
      final List<File> popPriorVcfFile = new CommandLineFiles(null, POPULATION_PRIORS, CommandLineFiles.EXISTS, CommandLineFiles.NOT_DIRECTORY, CommandLineFiles.TABIX).getFileList(mFlags);
      builder.populationPriors(popPriorVcfFile.get(0));
    }
    if (mFlags.isSet(MIN_VARIANT_ALLELIC_DEPTH)) {
      builder.minVariantAllelicDepth((Double) mFlags.getValue(MIN_VARIANT_ALLELIC_DEPTH));
    }
    if (mFlags.isSet(MIN_VARIANT_ALLELIC_FRACTION)) {
      builder.minVariantAllelicFraction((Double) mFlags.getValue(MIN_VARIANT_ALLELIC_FRACTION));
    }

    final AlleleBalanceFactor factor = mFlags.isSet(USE_ALLELIC_FRACTION)
      ? AlleleBalanceFactor.BINOMIAL
      : (AlleleBalanceFactor) mFlags.getValue(X_ALLELE_BALANCE_PROBABILITY);
    final AlleleBalanceProbability probability = getProbability(factor, (Double) mFlags.getValue(X_EXPECTED_ALLELIC_FRACTION));
    builder.alleleBalance(probability);

    builder.expandComplexReadQueries((Boolean) mFlags.getValue(X_EXPAND_COMPLEX_READ_QUERIES));
    builder.complexUseSoftClip((Boolean) mFlags.getValue(X_COMPLEX_USE_SOFT_CLIP));
    if (mFlags.isSet(X_FORCE_COMPLEX_REGION)) {
      builder.forceComplexRegion((String) mFlags.getValue(X_FORCE_COMPLEX_REGION));
    }

    // From here on pretty much needs the genome reader to be loaded
    final File genomeFile = (File) mFlags.getValue(CommonFlags.TEMPLATE_FLAG);
    SdfUtils.validateHasNames(genomeFile);
    final SequenceParams genomeParams = SequenceParams.builder().directory(genomeFile).mode(SequenceMode.UNIDIRECTIONAL).create();
    try {
      SdfUtils.validateNoDuplicates(genomeParams.reader(), false);
      builder.genome(genomeParams.readerParams());
      makeInputParams(builder, genomeParams, filterParams);
      return builder;
    } catch (final RuntimeException | IOException e) {
      genomeParams.close();
      throw e;
    }
  }

  private AlleleBalanceProbability getProbability(AlleleBalanceFactor factor, double expectedAlleleBalance) {
    switch (factor) {
      case BINOMIAL:
        return new BinomialAlleleBalance(expectedAlleleBalance);
      case HOEFFDING:
        return new HoeffdingAlleleBalance(expectedAlleleBalance);
      case NONE:
        return new NoAlleleBalance();
      default:
        throw new IllegalArgumentException();
    }
  }

  protected RegionRestriction getSimpleRegionRestriction() {
    if (mFlags.isSet(CommonFlags.RESTRICTION_FLAG)) {
      return new RegionRestriction((String) mFlags.getValue(CommonFlags.RESTRICTION_FLAG));
    } else {
      return null;
    }
  }

  private void makeInputParams(final VariantParamsBuilder builder, final SequenceParams genomeParams, SamFilterParams filterParams) throws IOException {
    final GenomeRelationships grf = grf();
    builder.genomeRelationships(grf);
    final Collection<File> inputFiles;
    if (mFlags.isSet(CommonFlags.RESTRICTION_FLAG) || mFlags.isSet(CommonFlags.BED_REGIONS_FLAG)) {
      inputFiles = new CommandLineFiles(CommonFlags.INPUT_LIST_FLAG, null, CommandLineFiles.EXISTS, CommandLineFiles.VARIANT_INPUT).getFileList(mFlags);
    } else {
      inputFiles = new CommandLineFiles(CommonFlags.INPUT_LIST_FLAG, null, CommandLineFiles.EXISTS).getFileList(mFlags);
    }
    final boolean useCalibration = !mFlags.isSet(MACHINE_ERRORS_FLAG) && !mFlags.isSet(NO_CALIBRATION);
    final SamCalibrationInputs inputs = new SamCalibrationInputs(inputFiles, useCalibration);
    final Collection<File> samFiles = inputs.getSamFiles();
    if (samFiles.isEmpty()) {
      throw new InvalidParamsException("No SAM files provided for input.");
    }
    Diagnostic.userLog("Input SAM files: " + samFiles);
    builder.mapped(samFiles);
    Calibrator c = null;
    if (useCalibration) {
      final Collection<File> calibrationFiles = inputs.getCalibrationFiles();
      Diagnostic.userLog("Input calibration files: " + calibrationFiles);
      builder.calibrations(calibrationFiles);
      if (calibrationFiles.isEmpty()) {
        throw new InvalidParamsException("No calibration files were found for input mappings. Please run 'rtg calibrate', or use --" + NO_CALIBRATION + " and set --" + COVERAGE_BYPASS_FLAG + " appropriately to proceed without calibration.");
      }
      if (calibrationFiles.size() != samFiles.size()) {
        throw new InvalidParamsException("Number of calibration files (" + calibrationFiles.size() + ") does not match number of SAM files (" + samFiles.size() + ").");
      }
      c = Calibrator.initCalibrator(calibrationFiles);
      if (c == null && mFlags.isSet(FILTER_DEPTH_MULTIPLIER_FLAG)) {
        throw new InvalidParamsException("Can not use --" + FILTER_DEPTH_MULTIPLIER_FLAG + " when calibration files are not being used.");
      }
      if (c == null && mFlags.isSet(COVERAGE_BYPASS_MULTIPLIER_FLAG)) {
        throw new InvalidParamsException("Can not use --" + COVERAGE_BYPASS_MULTIPLIER_FLAG + " when calibration files are not being used.");
      }
    }
    builder.ignoreIncompatibleSamHeaders(mFlags.isSet(X_IGNORE_SAM_HEADER_INCOMPATIBILITY));
    final SAMFileHeader uberHeader;
    if (mFlags.isSet(X_ALTERNATE_SAM_HEADER)) {
      final File altHeaderFile = (File) mFlags.getValue(X_ALTERNATE_SAM_HEADER);
      uberHeader = SamUtils.getSingleHeader(altHeaderFile);
    } else {
      uberHeader = SamUtils.getUberHeader(genomeParams.reader(), samFiles, mFlags.isSet(X_IGNORE_SAM_HEADER_INCOMPATIBILITY), grf == null ? null : grf.genomes());
    }
    builder.uberHeader(uberHeader);
    builder.referenceRanges(SamRangeUtils.createReferenceRanges(uberHeader, filterParams));

    if (c != null) {
      final RegionRestriction simpleRestriction = getSimpleRegionRestriction();
      final Map<String, String> readGroupToSampleId = SamUtils.getReadGroupToSampleId(uberHeader);
      final Map<String, Integer> sequenceLengthMap = c.hasLengths() ? c.getSequenceLengths() : Calibrator.getNonNSequenceLengthMap(genomeParams.reader(), simpleRestriction);
      final CalibratedPerSequenceExpectedCoverage expectedCoverages = new CalibratedPerSequenceExpectedCoverage(c, sequenceLengthMap, readGroupToSampleId, simpleRestriction);
      builder.expectedCoverage(expectedCoverages);
      if (!mFlags.isSet(FILTER_DEPTH_FLAG) && !mFlags.isSet(COVERAGE_BYPASS_FLAG)) {
        if (mFlags.isSet(FILTER_DEPTH_MULTIPLIER_FLAG)) {
          final ThresholdFunction filterFunc = mFlags.isSet(X_ALT_MULTIPLIER_FLAG) ? ThresholdFunction.SQRT_MULT : ThresholdFunction.SIMPLE_MULT; //ThresholdFunction.SIMPLE_MULT : ThresholdFunction.SQRT_MULT;
          final CalibratedPerSequenceThreshold filterThreshold = new CalibratedPerSequenceThreshold(expectedCoverages, (Double) mFlags.getValue(FILTER_DEPTH_MULTIPLIER_FLAG), filterFunc);
          builder.maxCoverageFilter(filterThreshold);
        }

        final ThresholdFunction processingFunc = mFlags.isSet(X_ALT_MULTIPLIER_FLAG) ? ThresholdFunction.SQRT_MULT : ThresholdFunction.SIMPLE_MULT; //ThresholdFunction.SIMPLE_MULT : ThresholdFunction.SQRT_MULT;
        final CalibratedPerSequenceThreshold processingThreshold = new CalibratedPerSequenceThreshold(expectedCoverages, (Double) mFlags.getValue(COVERAGE_BYPASS_MULTIPLIER_FLAG), processingFunc);
        builder.maxCoverageBypass(processingThreshold);
      }

      checkLowCoverage(expectedCoverages, new SexMemo(genomeParams.reader(), null));
    }
    if (c == null || mFlags.isSet(FILTER_DEPTH_FLAG) || mFlags.isSet(COVERAGE_BYPASS_FLAG)) {
      if (mFlags.isSet(FILTER_DEPTH_FLAG)) {
        builder.maxCoverageFilter(new StaticThreshold((Integer) mFlags.getValue(FILTER_DEPTH_FLAG)));
      }
      final Integer coverageBypassValue = (Integer) mFlags.getValue(COVERAGE_BYPASS_FLAG);
      final HashSet<String> samples = new HashSet<>();
      for (final SAMReadGroupRecord rgr : uberHeader.getReadGroups()) {
        samples.add(rgr.getSample());
      }
      builder.maxCoverageBypass(new StaticThreshold(coverageBypassValue, coverageBypassValue * Math.max(1, samples.size())));
    }
    builder.calibrator(c);
  }

  private void checkLowCoverage(CalibratedPerSequenceExpectedCoverage expectedCoverages, SexMemo sexMemo) {
    boolean doWarning = false;
    for (String sample : expectedCoverages.samples()) {
      for (String seq : expectedCoverages.sequences()) {
        if (sexMemo.isAutosome(seq)) { //only check sequences we expect normal coverage of
          final double coverage = expectedCoverages.expectedCoverage(seq, sample);
          if (coverage < LOW_COVERAGE_WARNING_THRESHOLD) {
            Diagnostic.warning(String.format("Calibration indicates very low coverage (%.1f) for sample: %s on sequence: %s", coverage, sample, seq));
            doWarning = true;
          }
        }
      }
    }
    if (doWarning) {
      Diagnostic.warning("Check that appropriate target regions were supplied during mapping/calibration");
    }
  }


  private void parseThreads(final VariantParamsBuilder builder) {
    builder.ioThreads(CommonFlags.parseThreads((Integer) mFlags.getValue(X_IO_THREADS)));
    builder.execThreads(CommonFlags.parseThreads((Integer) mFlags.getValue(CommonFlags.THREADS_FLAG)));
    if (mFlags.isSet(X_THREADING_ENVIRONMENT)) {
      final String val = (String) mFlags.getValue(X_THREADING_ENVIRONMENT);
      if ("single".equalsIgnoreCase(val)) {
        builder.threadingEnvironment(ThreadingEnvironment.SINGLE);
      } else if ("parallel".equalsIgnoreCase(val)) {
        builder.threadingEnvironment(ThreadingEnvironment.PARALLEL);
      } else {
        final String[] split = val.split("=");
        final Long l = Long.valueOf(split[1]);
        builder.threadingEnvironment(ThreadingEnvironment.RANDOM).threadingEnvironmentSeed(l);
      }
    }
  }

  private void processDefaultAndMaxMappingQuality(final VariantParamsBuilder builder) {
    final int qDefault = (Integer) mFlags.getValue(X_Q_DEFAULT_FLAG);
    if (qDefault < 0 || qDefault > 63) {
      throw new InvalidParamsException(ErrorType.INVALID_INTEGER, X_Q_DEFAULT_FLAG, "" + qDefault, "0", "63");
    }
    builder.defaultQuality(qDefault);
    final int rDefault = (Integer) mFlags.getValue(R_DEFAULT_FLAG);
    if (rDefault < 0 || rDefault > 63) {
      throw new InvalidParamsException(ErrorType.INVALID_INTEGER, R_DEFAULT_FLAG, "" + rDefault, "0", "63");
    }
    final int urDefault = (Integer) mFlags.getValue(UNMATED_R_DEFAULT_FLAG);
    if (urDefault < 0 || urDefault > 63) {
      throw new InvalidParamsException(ErrorType.INVALID_INTEGER, UNMATED_R_DEFAULT_FLAG, "" + urDefault, "0", "63");
    }
    builder.defaultMatedReadQuality(rDefault);
    builder.defaultUnmatedReadQuality(urDefault);
    final int rMax = (Integer) mFlags.getValue(X_R_MAX_FLAG);
    if (rMax < 0) {
      throw new InvalidParamsException(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, X_R_MAX_FLAG, "" + rMax, "0");
    }
    final int urMax = (Integer) mFlags.getValue(X_UNMATED_R_MAX_FLAG);
    if (urMax < 0) {
      throw new InvalidParamsException(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, X_UNMATED_R_MAX_FLAG, "" + urMax, "0");
    }
    builder.maxMatedReadQuality(rMax);
    builder.maxUnmatedReadQuality(urMax);
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  protected VariantStatistics getStatistics(VariantParams params) {
    return new VariantStatistics(params.directory());
  }

  @Override
  protected abstract IORunnable task(VariantParams params, OutputStream out) throws IOException;
}
