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

package com.rtg.variant.sv.discord.pattern;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ParamsCli;
import com.rtg.sam.SamFilterOptions;
import com.rtg.util.IORunnable;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.variant.sv.discord.DiscordantToolCli;

/**
 * Look for recognizable breakpoint patterns
 */
public class SvPatternsCli extends ParamsCli<BreakpointPatternParams> {

  /** Flag name for max fragment length */
  public static final String MAX_FRAGMENT_LENGTH = "max-fragment-length";
  /** Flag name for max same distance */
  public static final String MAX_SAME_DISTANCE = "max-same-distance";

  private static class PatternsValidator implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      return CommonFlags.validateOutputDirectory(flags)
      && CommonFlags.checkFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, Integer.MAX_VALUE)
      && CommonFlags.validateRegion(flags);
    }
  }

  @Override
  public String moduleName() {
    return "svpatterns";
  }

  @Override
  public String description() {
    return null;
  }

  @Override
  protected IORunnable task(BreakpointPatternParams params, OutputStream out) {
    return new SvPatternsTask(params, out);
  }

  @Override
  protected BreakpointPatternParams makeParams() throws InvalidParamsException, IOException {
    return makeParamsLocal(mFlags);
  }

  protected static BreakpointPatternParams makeParamsLocal(CFlags flags) throws IOException {
    final BreakpointPatternParams.Builder builder = BreakpointPatternParams.builder();
    if (flags.isSet(CommonFlags.RESTRICTION_FLAG)) {
      builder.region(new RegionRestriction((String) flags.getValue(CommonFlags.RESTRICTION_FLAG)));
    }
    if (flags.isSet(MAX_FRAGMENT_LENGTH)) {
      builder.fragmentLength((Integer) flags.getValue(MAX_FRAGMENT_LENGTH));
    }
    if (flags.isSet(MAX_SAME_DISTANCE)) {
      builder.sameDistance((Integer) flags.getValue(MAX_SAME_DISTANCE));
    }
    if (flags.isSet(DiscordantToolCli.MIN_SUPPORT_FLAG)) {
      builder.setMinDepth((Integer) flags.getValue(DiscordantToolCli.MIN_SUPPORT_FLAG));
    }
    return builder.directory((File) flags.getValue(CommonFlags.OUTPUT_FLAG))
        .files(CommonFlags.getFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, false))
        .create();
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  @Override
  protected void initFlags() {
    initLocalFlags(mFlags);
  }

  protected static void initLocalFlags(CFlags flags) {
    flags.setValidator(new PatternsValidator());
    CommonFlags.initOutputDirFlag(flags);
    CommonFlagCategories.setCategories(flags);
    SamFilterOptions.registerRestrictionFlag(flags);
    final Flag<File> inFlag = flags.registerRequired(File.class, CommonFlags.FILE, "VCF format files of discordant breakpoints");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    final Flag<File> listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, CommonFlags.FILE, "file containing a list of VCF format files (1 per line) of discordant breakpoints").setCategory(INPUT_OUTPUT);
    flags.addRequiredSet(inFlag);
    flags.addRequiredSet(listFlag);
    flags.registerOptional(MAX_FRAGMENT_LENGTH, Integer.class, CommonFlags.INT, "how far from the breakpoint to look ahead for inversions", BreakpointPatternParams.DEFAULT_FRAGMENT_LENGTH).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(MAX_SAME_DISTANCE, Integer.class, CommonFlags.INT, "how far apart can breakpoints be yet still be considered the same place", BreakpointPatternParams.DEFAULT_SAME_DISTANCE).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    DiscordantToolCli.registerMinSupport(flags);
  }

  /**
   * Command line
   * @param args command line args
   */
  public static void main(String[] args) {
    new SvPatternsCli().mainInit(args, System.out, System.err);
  }
}
