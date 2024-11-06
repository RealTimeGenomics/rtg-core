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
package com.rtg.variant.cnv.segment;

import static com.rtg.launcher.CommonFlags.FILE;
import static com.rtg.launcher.CommonFlags.FLOAT;
import static com.rtg.launcher.CommonFlags.INT;
import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.launcher.CommonFlags.STRING;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;
import static com.rtg.variant.cnv.segment.CnvPonBuildCli.NORMALIZED_COVERAGE_COLUMN;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.NoStatistics;
import com.rtg.launcher.OutputParams;
import com.rtg.launcher.ParamsCli;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.util.IORunnable;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;

/**
 * Provide segmentation entry point.
 */
public class SegmentCli extends ParamsCli<SegmentParams> {

  private static final String ALEPH_FLAG = "aleph";
  private static final String ALPHA_FLAG = "alpha";
  private static final String BETA_FLAG = "beta";
  private static final String MIN_SEGMENTS_FLAG = "Xmin-segments";
  private static final String MAX_SEGMENTS_FLAG = "Xmax-segments";

  private static final String MIN_LOGR_FLAG = "min-logr";
  private static final String MIN_BINS_FLAG = "min-bins";

  private static final String CASE_FLAG = "case";
  private static final String CONTROL_FLAG = "control";
  private static final String PANEL_FLAG = "panel";
  private static final String SAMPLE_FLAG = "sample";

  private static final String SUMMARY_FLAG = "summary-regions";

  private static final String COLUMN_FLAG = "Xcolumn";
  static final String GCBINS_FLAG = "Xgcbins";
  private static final String MIN_CASE_COV_FLAG = "min-case-coverage";
  private static final String MIN_CTRL_COV_FLAG = "min-control-coverage";
  private static final String MIN_NORM_CTRL_COV_FLAG = "min-norm-control-coverage";
  static final String COV_COLUMN_NAME = "coverage-column-name";
  static final String PANEL_COV_COLUMN_NAME = "Xpanel-coverage-column-name";
  private static final String GRAPHVIZ_SEGMENTATION = "Xgraphviz-segmentation";
  private static final String ABSORB_SINGLETONS_FLAG = "Xabsorb-singletons";

  static final String DEFAULT_COLUMN_NAME = "coverage";


  @Override
  public String moduleName() {
    return "segment";
  }

  @Override
  public String description() {
    return "segment depth of coverage data to identify copy number alterations";
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription("Segments depth of coverage data to identify regions of consistent copy number.");
    CommonFlagCategories.setCategories(mFlags);
    CommonFlags.initOutputDirFlag(mFlags);
    CommonFlags.initNoGzip(mFlags);
    CommonFlags.initIndexFlags(mFlags);
    CommonFlags.initReferenceTemplate(mFlags, true);
    mFlags.registerRequired(CASE_FLAG, File.class, FILE, "BED file supplying per-region coverage data for the sample").setCategory(INPUT_OUTPUT);
    final Flag<?> controlFlag = mFlags.registerOptional(CONTROL_FLAG, File.class, FILE, "BED file supplying per-region coverage data for control sample").setCategory(INPUT_OUTPUT);
    final Flag<?> panelFlag = mFlags.registerOptional(PANEL_FLAG, File.class, FILE, "BED file supplying per-region panel normalized coverage data").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(SUMMARY_FLAG, File.class, FILE, "BED file supplying gene-scale regions to report CNV interactions with").setCategory(INPUT_OUTPUT);

    mFlags.registerOptional(ALEPH_FLAG, Double.class, FLOAT, "weighting factor for inter-segment distances during energy scoring", 0.0).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(ALPHA_FLAG, Double.class, FLOAT, "weighting factor for intra-segment distances during energy scoring", 0.001).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(BETA_FLAG, Double.class, FLOAT, "segmentation sensitivity factor", 0.5).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(MIN_SEGMENTS_FLAG, Integer.class, INT, "lower bound on the number of segments to be produced", 1).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(MAX_SEGMENTS_FLAG, Integer.class, INT, "upper bound on the number of segments to be produced", Integer.MAX_VALUE).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(MIN_CASE_COV_FLAG, Double.class, FLOAT, "minimum case coverage required for a bin to be included in segmentation", 5.0).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(MIN_CTRL_COV_FLAG, Double.class, FLOAT, "minimum control coverage required for a bin to be included in segmentation").setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(MIN_NORM_CTRL_COV_FLAG, Double.class, FLOAT, "minimum normalized control or panel coverage required for a bin to be included in segmentation", 0.1).setCategory(SENSITIVITY_TUNING);

    mFlags.registerOptional(MIN_BINS_FLAG, Integer.class, INT, "minimum number of bins required for copy number alteration to be called", 1).setCategory(REPORTING);
    mFlags.registerOptional('r', MIN_LOGR_FLAG, Double.class, FLOAT, "minimum (absolute) log ratio required for copy number alteration to be called", 0.2).setCategory(REPORTING);
    mFlags.registerOptional('s', SAMPLE_FLAG, String.class, STRING, "sample name to associate with the case sample", "SAMPLE").setCategory(REPORTING);

    mFlags.registerOptional('c', COLUMN_FLAG, Integer.class, INT, "just run segmentation directly on the specified data column from the case file, ignoring control data", 0).setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(GCBINS_FLAG, Integer.class, INT, "number of bins when applying GC correction", 10).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(COV_COLUMN_NAME, String.class, STRING, "name of the coverage column in input data", DEFAULT_COLUMN_NAME).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(PANEL_COV_COLUMN_NAME, String.class, STRING, "name of the normalized coverage column in panel data", NORMALIZED_COVERAGE_COLUMN).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(GRAPHVIZ_SEGMENTATION, "if set, output a Graphviz file for viewing the segmentation tree for each chromosome").setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(ABSORB_SINGLETONS_FLAG, "absorb single bins into closest scoring adjacent segment").setCategory(SENSITIVITY_TUNING);
    mFlags.addRequiredSet(controlFlag);
    mFlags.addRequiredSet(panelFlag);
    mFlags.setValidator(flags -> CommonFlags.validateOutputDirectory(flags)
      && flags.checkInRange(MIN_LOGR_FLAG, 0, Double.MAX_VALUE)
      && flags.checkInRange(COLUMN_FLAG, 0, Integer.MAX_VALUE)
      && flags.checkInRange(MIN_CASE_COV_FLAG, 0, Double.MAX_VALUE)
      && flags.checkInRange(MIN_CTRL_COV_FLAG, 0, Double.MAX_VALUE)
      && flags.checkInRange(GCBINS_FLAG, 0, Integer.MAX_VALUE)
      && CommonFlags.validateInputFile(flags, CASE_FLAG)
      && CommonFlags.validateInputFile(flags, CONTROL_FLAG)
      && CommonFlags.validateInputFile(flags, SUMMARY_FLAG)
      && flags.checkInRange(MIN_SEGMENTS_FLAG, 1, Integer.MAX_VALUE)
      && flags.checkXor(COLUMN_FLAG, CONTROL_FLAG, PANEL_FLAG)
      && flags.checkNand(MIN_CTRL_COV_FLAG, MIN_NORM_CTRL_COV_FLAG)
      && checkMinMax(flags)
    );
  }

  private static boolean checkMinMax(final CFlags flags) {
    if ((Integer) flags.getValue(MIN_SEGMENTS_FLAG) > (Integer) flags.getValue(MAX_SEGMENTS_FLAG)) {
      flags.setParseMessage("--" + MIN_SEGMENTS_FLAG + " cannot exceed --" + MAX_SEGMENTS_FLAG);
      return false;
    }
    return true;
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(OUTPUT_FLAG);
  }

  @Override
  protected IORunnable task(SegmentParams params, OutputStream out) {
    return new SegmentTask(params, out, new NoStatistics());
  }

  @Override
  protected SegmentParams makeParams() throws IOException {
    final SegmentParams.SegmentParamsBuilder builder = SegmentParams.builder();
    final boolean gzip = !mFlags.isSet(CommonFlags.NO_GZIP);

    builder.name(mFlags.getName());
    final OutputParams outParams = new OutputParams(outputDirectory(), gzip);
    builder.outputParams(outParams);
    builder.genome(SequenceParams.builder().directory((File) mFlags.getValue(CommonFlags.TEMPLATE_FLAG)).mode(SequenceMode.UNIDIRECTIONAL).create().readerParams());

    if (mFlags.isSet(SUMMARY_FLAG)) {
      builder.summaryRegionsFile((File) mFlags.getValue(SUMMARY_FLAG));
    }
    if (mFlags.isSet(PANEL_FLAG)) {
      builder.panelFile((File) mFlags.getValue(PANEL_FLAG));
    }
    if (mFlags.isSet(CONTROL_FLAG)) {
      builder.controlFile((File) mFlags.getValue(CONTROL_FLAG));
    }
    builder.caseFile((File) mFlags.getValue(CASE_FLAG));

    if (mFlags.isSet(COLUMN_FLAG)) {
      builder.precomputedColumn((Integer) mFlags.getValue(COLUMN_FLAG));
    }
    builder.coverageColumnName((String) mFlags.getValue(COV_COLUMN_NAME));
    builder.panelCoverageColumnName((String) mFlags.getValue(PANEL_COV_COLUMN_NAME));
    if (mFlags.isSet(MIN_CTRL_COV_FLAG)) {
      builder.minControlCoverage((Double) mFlags.getValue(MIN_CTRL_COV_FLAG));
    }
    builder.minNormControlCoverage((Double) mFlags.getValue(MIN_NORM_CTRL_COV_FLAG));
    builder.sampleName((String) mFlags.getValue(SAMPLE_FLAG));
    builder.minCaseCoverage((Double) mFlags.getValue(MIN_CASE_COV_FLAG));
    builder.gcBins((Integer) mFlags.getValue(GCBINS_FLAG));
    builder.minBins((Integer) mFlags.getValue(MIN_BINS_FLAG));
    builder.absorbSingletons(mFlags.isSet(ABSORB_SINGLETONS_FLAG));
    builder.minLogR((Double) mFlags.getValue(MIN_LOGR_FLAG));
    builder.minSegments((Integer) mFlags.getValue(MIN_SEGMENTS_FLAG));
    builder.maxSegments((Integer) mFlags.getValue(MAX_SEGMENTS_FLAG));
    builder.aleph((Double) mFlags.getValue(ALEPH_FLAG));
    builder.alpha((Double) mFlags.getValue(ALPHA_FLAG));
    builder.beta((Double) mFlags.getValue(BETA_FLAG));
    builder.graphviz(mFlags.isSet(GRAPHVIZ_SEGMENTATION));
    return builder.create();
  }
}
