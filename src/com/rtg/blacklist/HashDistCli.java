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
package com.rtg.blacklist;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.launcher.BuildParams;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ParamsCli;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.ngs.MapFlags;
import com.rtg.usage.UsageMetric;
import com.rtg.util.IORunnable;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;

/**
 */
public class HashDistCli extends ParamsCli<HashDistParams> {

  private static final String MAX_COUNT_FLAG = "max-count";
  private static final String BLACKLIST_THRESHOLD_FLAG = "blacklist-threshold";
  private static final String INSTALL_BLACKLIST = "install-blacklist";
  private static final String HASHMAP_SIZE_FACTOR = "hashmap-size-factor";
  private static final String DESCRIPTION = "Counts the number of times hashes occur and produces a histogram. Optionally produces a blacklist of highly occuring hashes";

  @Override
  protected IORunnable task(HashDistParams params, OutputStream out) throws IOException {
    return new HashDist(params, out, new NullStatistics(), new UsageMetric());
  }

  @Override
  protected HashDistParams makeParams() throws InvalidParamsException, IOException {
    final BuildParams buildParams = BuildParams.builder()
      .windowSize((Integer) mFlags.getValue(MapFlags.WORDSIZE_FLAG))
      .stepSize((Integer) mFlags.getValue(MapFlags.STEP_FLAG))
      .sequences(
        SequenceParams.builder().directory((File) mFlags.getAnonymousValue(0)).mode(SequenceMode.UNIDIRECTIONAL).create()
      ).create();
    return HashDistParams.builder()
      .buildParams(buildParams)
      .installBlacklist(mFlags.isSet(INSTALL_BLACKLIST))
      .hashMapSizeFactor((Double) mFlags.getValue(HASHMAP_SIZE_FACTOR))
      .blacklistThreshold((Integer) mFlags.getValue(BLACKLIST_THRESHOLD_FLAG))
      .blacklist(mFlags.isSet(BLACKLIST_THRESHOLD_FLAG))
      .numberThreads(CommonFlags.parseThreads((Integer) mFlags.getValue(CommonFlags.THREADS_FLAG)))
      .directory((File) mFlags.getValue(CommonFlags.OUTPUT_FLAG))
      .threshold((Integer) mFlags.getValue(MAX_COUNT_FLAG))
      .create();
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  private static void initFlags(CFlags flags) {
    flags.registerExtendedHelp();
    flags.setDescription(DESCRIPTION);
    CommonFlagCategories.setCategories(flags);
    CommonFlags.initOutputDirFlag(flags);
    CommonFlags.initThreadsFlag(flags);
    flags.registerRequired(File.class, CommonFlags.SDF, "SDF data set to run against").setCategory(CommonFlagCategories.INPUT_OUTPUT).setMaxCount(1).setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerOptional('w', MapFlags.WORDSIZE_FLAG, Integer.class, CommonFlags.INT, "number of bases in hash", 22).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional('s', MapFlags.STEP_FLAG, Integer.class, CommonFlags.INT, "step size", 1).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(MAX_COUNT_FLAG, Integer.class, CommonFlags.INT, "soft minimum for hash count (i.e. will record exact counts of at least this much)", 500).setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(BLACKLIST_THRESHOLD_FLAG, Integer.class, CommonFlags.INT, "produces a blacklist of hashes with counts exceeding this value").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(INSTALL_BLACKLIST, "Adds blacklist to SDF for use in mapping").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerOptional(HASHMAP_SIZE_FACTOR, Double.class, CommonFlags.FLOAT, "Multiplier for the minimum size of the hashmap", 1.0).setCategory(CommonFlagCategories.UTILITY);
    flags.setValidator(new HashToolsCliValidator());
  }

  @Override
  public String moduleName() {
    return "hashdist";
  }

  @Override
  public String description() {
    return DESCRIPTION;
  }


  static final class HashToolsCliValidator implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      if (!CommonFlags.validateFlagBetweenValues(flags, MapFlags.WORDSIZE_FLAG, 1, 32)) {
        return false;
      }
      if (!CommonFlags.validateFlagBetweenValues(flags, HASHMAP_SIZE_FACTOR, 0.0, 5.0)) {
        return false;
      }
      if (!CommonFlags.validateSDF((File) flags.getAnonymousValue(0))) {
        return false;
      }
      if (!CommonFlags.validateThreads(flags)) {
        return false;
      }
      if (!CommonFlags.validateOutputDirectory(flags)) {
        return false;
      }
      if (flags.isSet(INSTALL_BLACKLIST) && !flags.isSet(BLACKLIST_THRESHOLD_FLAG)) {
        flags.setParseMessage("You must set --" + BLACKLIST_THRESHOLD_FLAG + " in order to install a blacklist");
        return false;
      }
      return true;
    }
  }


  /**
   * Main program for dusting. Use -h to get help.
   * @param args command line arguments.
   */
  public static void main(final String[] args) {
    new HashDistCli().mainExit(args);
  }
}
