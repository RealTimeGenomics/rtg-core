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
