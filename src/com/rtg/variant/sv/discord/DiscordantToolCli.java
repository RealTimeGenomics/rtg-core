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
package com.rtg.variant.sv.discord;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ParamsCli;
import com.rtg.launcher.SequenceParams;
import com.rtg.sam.SamFilterOptions;
import com.rtg.util.IORunnable;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.variant.sv.SvCliUtils;
import com.rtg.variant.sv.SvCliUtils.SvValidator;
import com.rtg.variant.sv.discord.DiscordantToolParams.DiscordantToolParamsBuilder;

/**
 */
public class DiscordantToolCli extends ParamsCli<DiscordantToolParams> {

  private static final String XDEBUG = "Xdebug-output";
  private static final String BED = "bed";
  private static final String MODULE_NAME = "discord";
  /** Flag name for minimum support */
  public static final String MIN_BREAKPOINT_DEPTH = "min-support";
  private static final String INTERSECTIONS = "consistent-only";
  /** Default value for <code>MIN_BREAKPOINT_DEPTH</code> */
  public static final Integer DEFAULT_MIN_DEPTH = 3;
  /** Description of min-support flag */
  public static final String MIN_SUPPORT_DESCRIPTION = "minimum number of supporting reads for a breakend";

  private static class DiscordantToolValidator extends SvValidator {

    @Override
    public boolean isValid(final CFlags flags) {
      if (!super.isValid(flags)) {
        return false;
      }
      if ((Integer) flags.getValue(MIN_BREAKPOINT_DEPTH) < 1) {
        Diagnostic.error(ErrorType.EXPECTED_POSITIVE, MIN_BREAKPOINT_DEPTH);
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
    initFlags(mFlags);
  }

  protected static void initFlags(CFlags flags) {
    SvCliUtils.initCommonFlags(flags);
    flags.setDescription("Analyses SAM records to determine the location of breakends.");
    flags.setValidator(new DiscordantToolValidator());
    SamFilterOptions.registerMaxHitsFlag(flags, 'c');
    flags.registerOptional('s', MIN_BREAKPOINT_DEPTH, Integer.class, "INT", MIN_SUPPORT_DESCRIPTION, DEFAULT_MIN_DEPTH).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(INTERSECTIONS, "only include breakends with internally consistent supporting reads").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(BED, "produce output in BED format in addition to VCF").setCategory(INPUT_OUTPUT);
//    flags.registerOptional(AbstractMultisampleCli.MAX_AMBIGUITY, IntegerOrPercentage.class, "int", "threshold for ambiguity above which calls are not made").setCategory(SENSITIVITY_TUNING);
//    flags.registerOptional(AbstractMultisampleCli.MAX_COVERAGE_FLAG, Integer.class, "int", "if set, will only output variants where coverage is less than this amount").setCategory(SENSITIVITY_TUNING);
    //X flags
    flags.registerOptional(XDEBUG, "produce debug output in addition to VCF").setCategory(INPUT_OUTPUT);
  }

  @Override
  protected DiscordantToolParams makeParams() throws InvalidParamsException, IOException {
    return makeParams(mFlags);
  }

  protected static DiscordantToolParams makeParams(CFlags flags) throws InvalidParamsException, IOException {
    final DiscordantToolParamsBuilder builder = DiscordantToolParams.builder();
    SvCliUtils.populateCommonParams(builder, SequenceParams.builder().useMemReader(false), flags);

    return builder.bedOutput(flags.isSet(BED))
        .minBreakpointDepth((Integer) flags.getValue(MIN_BREAKPOINT_DEPTH))
//        .maxCoverage((Integer) flags.getValue(AbstractMultisampleCli.MAX_COVERAGE_FLAG))
//        .maxAmbiguity(flags.isSet(AbstractMultisampleCli.MAX_AMBIGUITY) ? (((IntegerOrPercentage) flags.getValue(AbstractMultisampleCli.MAX_AMBIGUITY)).getValue(100) / 100.0) : null)
        .intersectionOnly(flags.isSet(INTERSECTIONS))
        .debugOutput(flags.isSet(XDEBUG))
        .create();
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  @Override
  protected IORunnable task(final DiscordantToolParams params, final OutputStream out) throws IOException {
    return new DiscordantTool(params, out);
  }
  /**
   * @param args command line arguments
   */
  public static void main(String[] args) {
    new DiscordantToolCli().mainExit(args);
  }

}
