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
package com.rtg.blacklist;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.launcher.BuildParams;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.NoStatistics;
import com.rtg.launcher.ParamsCli;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.mode.SequenceType;
import com.rtg.ngs.MapFlags;
import com.rtg.reader.ReaderUtils;
import com.rtg.usage.UsageMetric;
import com.rtg.util.IORunnable;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;

/**
 * Command-line entry point for the <code>hashdist</code> command.
 */
public class HashDistCli extends ParamsCli<HashDistParams> {

  private static final String MAX_COUNT_FLAG = "max-count";
  private static final String BLACKLIST_THRESHOLD_FLAG = "blacklist-threshold";
  private static final String INSTALL_BLACKLIST = "install-blacklist";
  private static final String HASHMAP_SIZE_FACTOR = "hashmap-size-factor";

  @Override
  protected IORunnable task(HashDistParams params, OutputStream out) {
    return new HashDist(params, out, new NoStatistics(), new UsageMetric());
  }

  @Override
  protected HashDistParams makeParams() throws IOException {
    final File sdfDir = (File) mFlags.getAnonymousValue(0);
    if (ReaderUtils.isPairedEndDirectory(sdfDir)) {
      throw new InvalidParamsException("Analysis of paired-end SDFs is not supported");
    }
    final SequenceParams sequences = SequenceParams.builder().directory(sdfDir).mode(SequenceMode.UNIDIRECTIONAL).create();
    if (sequences.reader().type() != SequenceType.DNA) {
      throw new InvalidParamsException("Analysis of non-DNA containing SDFs is not supported");
    }
    final BuildParams buildParams = BuildParams.builder()
      .windowSize((Integer) mFlags.getValue(MapFlags.WORDSIZE_FLAG))
      .stepSize((Integer) mFlags.getValue(MapFlags.STEP_FLAG))
      .sequences(sequences).create();
    return HashDistParams.builder()
      .buildParams(buildParams)
      .installBlacklist(mFlags.isSet(INSTALL_BLACKLIST))
      .hashMapSizeFactor((Double) mFlags.getValue(HASHMAP_SIZE_FACTOR))
      .blacklistThreshold((Integer) mFlags.getValue(BLACKLIST_THRESHOLD_FLAG))
      .makeBlacklist(mFlags.isSet(BLACKLIST_THRESHOLD_FLAG))
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
    flags.setDescription("Counts the number of times k-mers occur in an SDF and produces a histogram. Optionally creates a blacklist of highly occurring hashes that can be used to increase mapping speed.");
    CommonFlagCategories.setCategories(flags);
    CommonFlags.initOutputDirFlag(flags);
    CommonFlags.initThreadsFlag(flags);
    flags.registerRequired(File.class, CommonFlags.SDF, "SDF containing sequence data to analyse").setCategory(CommonFlagCategories.INPUT_OUTPUT).setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerOptional('w', MapFlags.WORDSIZE_FLAG, Integer.class, CommonFlags.INT, "number of bases in each hash", 22).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional('s', MapFlags.STEP_FLAG, Integer.class, CommonFlags.INT, "step size", 1).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(MAX_COUNT_FLAG, Integer.class, CommonFlags.INT, "soft minimum for hash count (i.e. will record exact counts of at least this value)", 500).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(BLACKLIST_THRESHOLD_FLAG, Integer.class, CommonFlags.INT, "if set, output a blacklist containing all k-mer hashes with counts exceeding this value").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(HASHMAP_SIZE_FACTOR, Double.class, CommonFlags.FLOAT, "multiplier for the minimum size of the hash map", 1.0).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(INSTALL_BLACKLIST, "install the blacklist into the SDF for use during mapping").setCategory(CommonFlagCategories.UTILITY);
    flags.setValidator(new HashToolsCliValidator());
  }

  @Override
  public String moduleName() {
    return "hashdist";
  }

  @Override
  public String description() {
    return "analyse the k-mer distribution within an SDF";
  }


  static final class HashToolsCliValidator implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      if (!flags.checkInRange(MapFlags.WORDSIZE_FLAG, 1, 32)) {
        return false;
      }
      if (!flags.checkInRange(HASHMAP_SIZE_FACTOR, 0.0, 5.0)) {
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

}
