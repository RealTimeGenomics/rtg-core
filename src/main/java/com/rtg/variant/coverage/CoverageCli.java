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

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.OutputParams;
import com.rtg.launcher.ParamsCli;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.sam.SamFilterOptions;
import com.rtg.util.IORunnable;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.diagnostic.Diagnostic;

/**
 */
public class CoverageCli extends ParamsCli<CoverageParams> {

  private static final String SMOOTHING_LEVEL_FLAG = "smoothing";
  private static final String ERROR_RATES_FLAG = "Xerror-rates";
  private static final String BEDGRAPH_FLAG = "bedgraph";
  private static final String PER_BASE_FLAG = "per-base";
  private static final String PER_REGION_FLAG = "per-region";
  private static final String X_COVERAGE_THRESHOLD_FLAG = "Xcoverage-threshold";
  private static final String X_FOLD_PCT_FLAG = "Xfold-penalty-target-percent";
  private static final String X_IGNORE_SAM_HEADER_INCOMPATIBILITY_FLAG = "Xignore-incompatible-sam-headers";
  private static final String X_BINARIZE_BED_FLAG = "Xbinarize-bed";
  private static final String X_CALLABILITY_FLAG = "Xcallability";
  private static final String X_DISABLE_HTML_REPORT_FLAG = "Xdisable-html-report";

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
    mFlags.setDescription("Measures and reports coverage depth of read alignments across a reference.");
    CommonFlagCategories.setCategories(mFlags);
    final Flag<File> inFlag = mFlags.registerRequired(File.class, CommonFlags.FILE, "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    final Flag<File> listFlag = mFlags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, CommonFlags.FILE, "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);
    CommonFlags.initOutputDirFlag(mFlags);
    CommonFlags.initReferenceTemplate(mFlags, false);
    CommonFlags.initNoGzip(mFlags);
    mFlags.registerOptional(PER_BASE_FLAG, "if set, output per-base counts in TSV format (suppresses BED file output)").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(PER_REGION_FLAG, "if set, output BED/BEDGRAPH entries per-region rather than every coverage level change").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(BEDGRAPH_FLAG, "if set, output in BEDGRAPH format (suppresses BED file output)").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional('s', SMOOTHING_LEVEL_FLAG, Integer.class, CommonFlags.INT, "smooth with this number of neighboring values (0 means no smoothing)", 50).setCategory(SENSITIVITY_TUNING);
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
    mFlags.registerOptional(X_COVERAGE_THRESHOLD_FLAG, Integer.class, CommonFlags.INT, "coverage threshold for breadth computation (and binarization, if enabled)", 1).setCategory(REPORTING);
    mFlags.registerOptional(X_BINARIZE_BED_FLAG, "if set, binarize BED outputs").setCategory(REPORTING);
    mFlags.registerOptional(X_CALLABILITY_FLAG, Integer.class, CommonFlags.INT, "report callability with respect to this minimum coverage level").setCategory(REPORTING);
    mFlags.registerOptional(X_FOLD_PCT_FLAG, Integer.class, CommonFlags.INT, "target percent when computing fold penalty", 80).setCategory(REPORTING);
    mFlags.registerOptional(X_IGNORE_SAM_HEADER_INCOMPATIBILITY_FLAG, "ignore incompatible SAM headers when merging SAM results").setCategory(UTILITY);
    mFlags.registerOptional(X_DISABLE_HTML_REPORT_FLAG, "disable HTML report output").setCategory(REPORTING);
    CommonFlags.initIndexFlags(mFlags);
    mFlags.addRequiredSet(inFlag);
    mFlags.addRequiredSet(listFlag);

    mFlags.setValidator(flags -> CommonFlags.validateOutputDirectory(flags)
      && CommonFlags.validateTemplate(flags)
      /*
       * The 8192 upper limit is based on the chunk size (10,000 as set in CoverageParams).
       * This value must be less than the chunk size since the window extends
       * this value either side of the point being smoothed, but must not cross two chunk boundaries.
       */
      && flags.checkInRange(SMOOTHING_LEVEL_FLAG, 0, 8192)
      && CommonFlags.checkFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, Integer.MAX_VALUE)
      && CommonFlags.validateThreads(flags)
      && SamFilterOptions.validateFilterFlags(flags, false)
      && flags.checkNand(PER_BASE_FLAG, BEDGRAPH_FLAG)
      && flags.checkNand(PER_BASE_FLAG, PER_REGION_FLAG)
      && flags.checkNand(PER_BASE_FLAG, SMOOTHING_LEVEL_FLAG)
      && flags.checkNand(PER_REGION_FLAG, X_BINARIZE_BED_FLAG)
      && flags.checkNand(PER_REGION_FLAG, X_CALLABILITY_FLAG)
      && flags.checkNand(PER_REGION_FLAG, SMOOTHING_LEVEL_FLAG)
    );
  }

  @Override
  protected CoverageParams makeParams() throws IOException {
    final CoverageParams.CoverageParamsBuilder builder = CoverageParams.builder();
    builder.name(mFlags.getName());
    final OutputParams outParams = new OutputParams((File) mFlags.getValue(CommonFlags.OUTPUT_FLAG), !mFlags.isSet(NO_GZIP));
    builder.outputParams(outParams);
    if (mFlags.isSet(CommonFlags.TEMPLATE_FLAG)) {
      builder.genome(SequenceParams.builder().directory((File) mFlags.getValue(CommonFlags.TEMPLATE_FLAG)).mode(SequenceMode.UNIDIRECTIONAL).create().readerParams());
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
    builder.foldTargetPercent((Integer) mFlags.getValue(X_FOLD_PCT_FLAG));
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
    return new com.rtg.variant.coverage.CoverageTask(params, out, new CoverageStatistics(outputDirectory(), params.disableHtmlReport(), params.foldTargetPercent()));
  }
}
