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

import static com.rtg.launcher.CommonFlags.FLOAT;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ParamsCli;
import com.rtg.launcher.SequenceParams;
import com.rtg.sam.SamFilterOptions;
import com.rtg.util.IORunnable;
import com.rtg.util.cli.CFlags;
import com.rtg.variant.sv.SvCliUtils;
import com.rtg.variant.sv.SvCliUtils.SvValidator;
import com.rtg.variant.sv.discord.DiscordantToolParams.DiscordantToolParamsBuilder;

/**
 */
public class DiscordantToolCli extends ParamsCli<DiscordantToolParams> {

  /** Flag name for minimum support */
  public static final String MIN_SUPPORT_FLAG = "min-support";
  private static final String INTERSECTIONS_FLAG = "consistent-only";
  private static final String BED_FLAG = "bed";
  private static final String OVERLAP_FRACTION = "overlap-fraction";
  private static final String NUM_DEVIATIONS = "Xnum-deviations"; // We may want to use an alternate model of gap distribution
  private static final String BAM_FLAG = "Xbam-output";
  private static final String DEBUG_FLAG = "Xdebug-output";
  private static final String MULTISAMPLE_FLAG = "Xmultisample";


  @Override
  public String moduleName() {
    return "discord";
  }

  @Override
  public String description() {
    return "detect structural variant breakends using discordant reads";
  }

  @Override
  protected void initFlags() {
    SvCliUtils.initCommonFlags(mFlags);
    mFlags.setDescription("Analyses SAM records to determine the location of breakends.");
    SamFilterOptions.registerMaxHitsFlag(mFlags, 'c');
    registerMinSupport(mFlags);
    mFlags.registerOptional(INTERSECTIONS_FLAG, "only include breakends with internally consistent supporting reads").setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(OVERLAP_FRACTION, Double.class, FLOAT, "assume this fraction of an aligned ready may may overlap a breakend", 0.01).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(NUM_DEVIATIONS, Double.class, FLOAT, "the number of standard deviations either side which are considered concordant", 4.0).setCategory(SENSITIVITY_TUNING);
    mFlags.registerOptional(BED_FLAG, "produce output in BED format in addition to VCF").setCategory(INPUT_OUTPUT);

    mFlags.registerOptional(BAM_FLAG, DiscordantTool.BamType.class, CommonFlags.STRING, "produce BAM output containing clusters of discordant alignments", DiscordantTool.BamType.NONE).setCategory(REPORTING);
    mFlags.registerOptional(DEBUG_FLAG, "produce debug output in addition to VCF").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(MULTISAMPLE_FLAG, "allow running on pooled data from multiple samples").setCategory(INPUT_OUTPUT);
    mFlags.setValidator(new SvValidator() {
      @Override
      public boolean isValid(final CFlags flags) {
        return super.isValid(flags)
          && flags.checkInRange(OVERLAP_FRACTION, 0.0, 1.0)
          && flags.checkInRange(NUM_DEVIATIONS, 3.0, Integer.MAX_VALUE)
          && flags.checkInRange(MIN_SUPPORT_FLAG, 1, Integer.MAX_VALUE);
      }
    });
  }

  /**
   * Registers a flag for the minimum number of reads that must support the breakend
   * @param flags flags to register with
   */
  public static void registerMinSupport(CFlags flags) {
    flags.registerOptional('s', MIN_SUPPORT_FLAG, Integer.class, CommonFlags.INT, "minimum number of supporting reads for a breakend", 3).setCategory(SENSITIVITY_TUNING);
  }

  @Override
  protected DiscordantToolParams makeParams() throws IOException {
    final DiscordantToolParamsBuilder builder = DiscordantToolParams.builder();
    SvCliUtils.populateCommonParams(builder, SequenceParams.builder(), mFlags);
    return builder.bedOutput(mFlags.isSet(BED_FLAG))
        .minBreakpointDepth((Integer) mFlags.getValue(MIN_SUPPORT_FLAG))
        .overlapFraction((Double) mFlags.getValue(OVERLAP_FRACTION))
        .numDeviations((Double) mFlags.getValue(NUM_DEVIATIONS))
        .intersectionOnly(mFlags.isSet(INTERSECTIONS_FLAG))
        .debugOutput(mFlags.isSet(DEBUG_FLAG))
        .allowMultisample(mFlags.isSet(MULTISAMPLE_FLAG))
        .bamOutput((DiscordantTool.BamType) mFlags.getValue(BAM_FLAG))
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
}
