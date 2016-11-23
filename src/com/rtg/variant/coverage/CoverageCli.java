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
package com.rtg.variant.coverage;

import static com.rtg.launcher.CommonFlags.NO_GZIP;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collection;

import com.rtg.launcher.BuildCommon;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.OutputParams;
import com.rtg.launcher.ParamsCli;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.sam.SamFilterOptions;
import com.rtg.util.IORunnable;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;

/**
 */
public class CoverageCli extends ParamsCli<CoverageParams> {

  private static final String TEMPLATE_FLAG = "template";
  private static final String SMOOTHING_LEVEL_FLAG = "smoothing";
  private static final String ERROR_RATES_FLAG = "Xerror-rates";
  private static final String BEDGRAPH_FLAG = "bedgraph";
  private static final String PER_BASE_FLAG = "per-base";
  private static final String PER_REGION_FLAG = "per-region";
  private static final String X_COVERAGE_THRESHOLD_FLAG = "Xcoverage-threshold";
  private static final String X_IGNORE_SAM_HEADER_INCOMPATIBILITY_FLAG = "Xignore-incompatible-sam-headers";
  private static final String X_BINARIZE_BED_FLAG = "Xbinarize-bed";
  private static final String X_CALLABILITY_FLAG = "Xcallability";
  private static final String X_DISABLE_HTML_REPORT_FLAG = "Xdisable-html-report";

  private static class CoverageValidator implements Validator {

    @Override
    public boolean isValid(final CFlags flags) {
      if (!CommonFlags.validateOutputDirectory(flags)) {
        return false;
      }

      if (flags.isSet(TEMPLATE_FLAG)) {
        final String genomePath = flags.getFlag(TEMPLATE_FLAG).getValue().toString();
        final File genome = new File(genomePath);
        if (!genome.exists()) {
          Diagnostic.error(ErrorType.INFO_ERROR, "The specified SDF, \"" + genome.getPath() + "\", does not exist.");
          return false;
        }
        if (!genome.isDirectory()) {
          Diagnostic.error(ErrorType.INFO_ERROR, "The specified file, \"" + genome.getPath() + "\", is not an SDF.");
          return false;
        }
      }

      if (flags.isSet(SMOOTHING_LEVEL_FLAG)) {
        /*
         * The 8192 upper limit is based on the chunk size (10,000 as set in CoverageParams).
         * This value must be less than the chunk size since the window extends
         * this value either side of the point being smoothed, but must not cross two chunk boundaries.
         */
        if (!flags.checkInRange(SMOOTHING_LEVEL_FLAG, 0, 8192)) {
          return false;
        }
      }

      if (!CommonFlags.checkFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, Integer.MAX_VALUE)) {
        return false;
      }

      if (!CommonFlags.validateThreads(flags)) {
        return false;
      }

      if (!SamFilterOptions.validateFilterFlags(flags, false)) {
        return false;
      }

      if (!flags.checkNand(PER_BASE_FLAG, BEDGRAPH_FLAG)) {
        return false;
      }
      if (!flags.checkNand(PER_BASE_FLAG, PER_REGION_FLAG)) {
        return false;
      }
      if (!flags.checkNand(PER_BASE_FLAG, SMOOTHING_LEVEL_FLAG)) {
        return false;
      }
      if (!flags.checkNand(PER_REGION_FLAG, X_BINARIZE_BED_FLAG)) {
        return false;
      }
      if (!flags.checkNand(PER_REGION_FLAG, X_CALLABILITY_FLAG)) {
        return false;
      }
      if (!flags.checkNand(PER_REGION_FLAG, SMOOTHING_LEVEL_FLAG)) {
        return false;
      }

      return true;
    }

  }

  @Override
  public String moduleName() {
    return "coverage";
  }

  @Override
  public String description() {
    return "calculate depth of coverage from SAM/BAM files";
  }

  @Override
  protected void initFlags() {
    mFlags.registerExtendedHelp();
    mFlags.setDescription("Measures and reports coverage depth of read alignments across a reference.");
    CommonFlagCategories.setCategories(mFlags);
    mFlags.setValidator(new CoverageValidator());
    final Flag inFlag = mFlags.registerRequired(File.class, "FILE", "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    final Flag listFlag = mFlags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, "FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);
    CommonFlags.initOutputDirFlag(mFlags);
    mFlags.registerOptional('t', TEMPLATE_FLAG, File.class, "SDF", "SDF of the reference genome the reads have been mapped against").setCategory(INPUT_OUTPUT);
    CommonFlags.initNoGzip(mFlags);
    mFlags.registerOptional(PER_BASE_FLAG, "if set, output per-base counts in TSV format (suppresses BED file output)").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(PER_REGION_FLAG, "if set, output BED/BEDGRAPH entries per-region rather than every coverage level change").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(BEDGRAPH_FLAG, "if set, output in BEDGRAPH format (suppresses BED file output)").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional('s', SMOOTHING_LEVEL_FLAG, Integer.class, "INT", "smooth with this number of neighboring values (0 means no smoothing)", 50).setCategory(SENSITIVITY_TUNING);
    CommonFlags.initThreadsFlag(mFlags);
    SamFilterOptions.registerMinMapQFlag(mFlags);
    SamFilterOptions.registerMaxHitsFlag(mFlags, 'c');
    SamFilterOptions.registerMaxASMatedFlag(mFlags, 'm');
    SamFilterOptions.registerMaxASUnmatedFlag(mFlags, 'u');
    SamFilterOptions.registerExcludeMatedFlag(mFlags);
    SamFilterOptions.registerExcludeUnmatedFlag(mFlags);
    SamFilterOptions.registerRestrictionFlag(mFlags);
    SamFilterOptions.registerBedRestrictionFlag(mFlags);
    SamFilterOptions.registerKeepDuplicatesFlag(mFlags);
    mFlags.registerOptional(ERROR_RATES_FLAG, "report statistics about sequencer error rates").setCategory(REPORTING);
    mFlags.registerOptional(X_COVERAGE_THRESHOLD_FLAG, Integer.class, "INT", "coverage threshold for breadth computation (and binarization, if enabled)", 1).setCategory(REPORTING);
    mFlags.registerOptional(X_BINARIZE_BED_FLAG, "if set, binarize BED outputs").setCategory(REPORTING);
    mFlags.registerOptional(X_CALLABILITY_FLAG, Integer.class, "INT", "report callability with respect to this minimum coverage level").setCategory(REPORTING);
    mFlags.registerOptional(X_IGNORE_SAM_HEADER_INCOMPATIBILITY_FLAG, "ignore incompatible SAM headers when merging SAM results").setCategory(UTILITY);
    mFlags.registerOptional(X_DISABLE_HTML_REPORT_FLAG, "disable HTML report output").setCategory(REPORTING);
    CommonFlags.initIndexFlags(mFlags);
    mFlags.addRequiredSet(inFlag);
    mFlags.addRequiredSet(listFlag);
  }

  @Override
  protected CoverageParams makeParams() throws InvalidParamsException, IOException {
    final CoverageParams.CoverageParamsBuilder builder = CoverageParams.builder();
    builder.name(mFlags.getName());
    final OutputParams outParams = new OutputParams((File) mFlags.getValue(CommonFlags.OUTPUT_FLAG), mFlags.isSet(BuildCommon.PROGRESS_FLAG), !mFlags.isSet(NO_GZIP));
    builder.outputParams(outParams);
    if (mFlags.isSet(TEMPLATE_FLAG)) {
      builder.genome(SequenceParams.builder().directory((File) mFlags.getValue(TEMPLATE_FLAG)).mode(SequenceMode.UNIDIRECTIONAL).create().readerParams());
    }
    final Collection<File> inputFiles = CommonFlags.getFileList(mFlags, CommonFlags.INPUT_LIST_FLAG, null, false);
    Diagnostic.userLog("Input SAM files: " + inputFiles);
    builder.mapped(inputFiles);
    builder.smoothing((Integer) mFlags.getValue(SMOOTHING_LEVEL_FLAG));
    builder.errorRates(mFlags.isSet(ERROR_RATES_FLAG));
    builder.tsvOutput(mFlags.isSet(PER_BASE_FLAG));
    builder.perRegion(mFlags.isSet(PER_REGION_FLAG));
    builder.bedgraphOutput(mFlags.isSet(BEDGRAPH_FLAG));
    builder.ioThreads(CommonFlags.parseIOThreads((Integer) mFlags.getValue(CommonFlags.THREADS_FLAG)));
    builder.outputIndex(!mFlags.isSet(CommonFlags.NO_INDEX));
    builder.ignoreIncompatibleSamHeaders(mFlags.isSet(X_IGNORE_SAM_HEADER_INCOMPATIBILITY_FLAG));
    builder.minimumCoverageThreshold((Integer) mFlags.getValue(X_COVERAGE_THRESHOLD_FLAG));
    builder.binarizeBed(mFlags.isSet(X_BINARIZE_BED_FLAG));
    if (mFlags.isSet(X_CALLABILITY_FLAG)) {
      builder.includeDeletions(true).binarizeBed(true).minimumCoverageThreshold((Integer) mFlags.getValue(X_CALLABILITY_FLAG));
    }
    if (mFlags.isSet(X_DISABLE_HTML_REPORT_FLAG)) {
      builder.disableHtmlReport(true);
    }
    return builder.filterParams(SamFilterOptions.makeFilterParamsBuilder(mFlags).excludeUnmapped(true).excludeUnplaced(true).excludeVariantInvalid(true).create()).create();
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  @Override
  protected IORunnable task(final CoverageParams params, final OutputStream out) {
    return new com.rtg.variant.coverage.CoverageTask(params, out, new CoverageStatistics(outputDirectory(), params.disableHtmlReport()));
  }
}
