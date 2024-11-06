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

package com.rtg.variant.sv;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.OutputParams;
import com.rtg.launcher.SequenceParams.SequenceParamsBuilder;
import com.rtg.sam.SamBamConstants;
import com.rtg.sam.SamFilterOptions;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.Utils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.sv.SvParams.SvParamsBuilder;

/**
 * Utility class for common structural variant command line interface code
 */
public final class SvCliUtils {

  private static final String RG_STATS_LIST_FLAG = "readgroup-stats-list-file";
  static final String RG_STATS_FILE = "readgroup-stats";
  private static final String RELABEL_FLAG = "readgroup-labels";

  private SvCliUtils() { }

  /**
   * <code>Validator</code> for common structural variant flags
   */
  public abstract static class SvValidator implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      return CommonFlags.validateOutputDirectory(flags)
        && CommonFlags.validateSDF(flags, CommonFlags.TEMPLATE_FLAG)
        && CommonFlags.checkFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, Integer.MAX_VALUE)
        && CommonFlags.checkFileList(flags, RG_STATS_LIST_FLAG, RG_STATS_FILE, Integer.MAX_VALUE)
        && CommonFlags.validateThreads(flags)
        && SamFilterOptions.validateFilterFlags(flags, false);
    }
  }

  /**
   * Initialise the flags common to all structural variant commands
   * @param flags the flags to initialise
   */
  public static void initCommonFlags(CFlags flags) {
    CommonFlagCategories.setCategories(flags);
    CommonFlags.initReferenceTemplate(flags, true);
    final Flag<File> inFlag = flags.registerRequired(File.class, CommonFlags.FILE, "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    final Flag<File> listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, CommonFlags.FILE, "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);
    CommonFlags.initOutputDirFlag(flags);
    CommonFlags.initNoGzip(flags);
    final Flag<File> rgstats = flags.registerOptional('r', RG_STATS_FILE, File.class, CommonFlags.FILE, "text file containing read group stats").setCategory(INPUT_OUTPUT);
    rgstats.setMinCount(0);
    rgstats.setMaxCount(Integer.MAX_VALUE);
    final Flag<File> rgListFlag = flags.registerOptional('R', RG_STATS_LIST_FLAG, File.class, CommonFlags.FILE, "file containing list of read group stats files (1 per line)").setCategory(INPUT_OUTPUT);
    CommonFlags.initThreadsFlag(flags);
    SamFilterOptions.registerMaxASMatedFlag(flags, 'm');
    SamFilterOptions.registerMaxASUnmatedFlag(flags, 'u');
    SamFilterOptions.registerRestrictionFlag(flags);
    flags.addRequiredSet(inFlag, rgstats);
    flags.addRequiredSet(listFlag, rgListFlag);
  }

  /**
   * Initialise the read-group label flag
   * @param flags the flags to initialise
   */
  public static void initRelabelFlag(CFlags flags) {
    flags.registerOptional(RELABEL_FLAG, File.class, CommonFlags.FILE, "file containing read group relabel mappings (1 per line)").setCategory(INPUT_OUTPUT);
  }

  /**
   * Populate a builder with the common flag settings
   * @param builder the builder to populate
   * @param genomeBuilder <code>SequenceParamsBuilder</code> to add the directory to and use
   * @param flags the flags to use to populate the builder
   * @throws InvalidParamsException if there are errors in the values of the command line flags
   * @throws IOException If an I/O error occurs
   */
  public static void populateCommonParams(SvParamsBuilder<?> builder, SequenceParamsBuilder genomeBuilder, CFlags flags) throws IOException {
    builder.name(flags.getName())
           .genome(genomeBuilder.directory((File) flags.getValue(CommonFlags.TEMPLATE_FLAG)).create().readerParams())
           .outputParams(new OutputParams((File) flags.getValue(CommonFlags.OUTPUT_FLAG), !flags.isSet(CommonFlags.NO_GZIP)))
           .ioThreads(CommonFlags.parseIOThreads((Integer) flags.getValue(CommonFlags.THREADS_FLAG)))
           .mapped(CommonFlags.getFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, false))
           .filterParams(SamFilterOptions.makeFilterParamsBuilder(flags)
             .requireSetFlags(SamBamConstants.SAM_READ_IS_PAIRED)
             .requireUnsetFlags(SamBamConstants.SAM_READ_IS_UNMAPPED | SamBamConstants.SAM_SECONDARY_ALIGNMENT | SamBamConstants.SAM_SUPPLEMENTARY_ALIGNMENT)
             .excludeDuplicates(true).excludeUnplaced(true).create());

    populateReadGroupStats(builder, flags);
  }

  static void populateReadGroupStats(SvParamsBuilder<?> builder, CFlags flags) throws IOException {
    final Map<String, String> relabelMap;
    if (flags.isSet(RELABEL_FLAG)) {
      relabelMap = ReadGroupStats.loadRelabelFile((File) flags.getValue(RELABEL_FLAG));
    } else {
      relabelMap = null;
    }
    builder.readGroupLabels(relabelMap);
    final Collection<File> statsFiles = CommonFlags.getFileList(flags, RG_STATS_LIST_FLAG, RG_STATS_FILE, false);
    final Map<String, ReadGroupStats> readGroupStats = ReadGroupStats.loadReadGroupStats(relabelMap, statsFiles.toArray(new File[0]));
    for (final ReadGroupStats e : readGroupStats.values()) {
      if (!e.isValid()) {
        throw new NoTalkbackSlimException("Invalid read group statistics for read group: " + e.id());
      //the below checks are what was previously being asserted upon output of rg stats file.
      } else if (e.meanLength() <= 10.0) {
        throw new NoTalkbackSlimException("Mean length in statistics for read group: " + e.id() + " must be greater than 10.0");
      } else if (e.fragmentMean() <= 20.0) {
        throw new NoTalkbackSlimException("Fragment mean in statistics for read group: " + e.id() + " must be greater than 20.0");
      } else if (e.fragmentStdDev() == 0.0 || e.gapStdDev() == 0.0 || e.properRate() == 0.0 || e.properRandomRate() == 0.0 || e.discordantRate() == 0.0 || e.unmatedRate() == 0.0) {
        throw new NoTalkbackSlimException("Statistics for read group: " + e.id() + " must be greater than 0.0");
      }
      warnStatistic(e.id(), "mean read length", e.meanLength());
      warnStatistic(e.id(), "fragment length mean", e.fragmentMean());
      warnStatistic(e.id(), "gap mean", e.gapMean());
    }
    builder.readGroupStatistics(readGroupStats);
  }

  private static void warnStatistic(String id, String msg, double x) {
    if (x > 10000.0) {
      Diagnostic.warning("Read group " + id + " statistics " + msg + " (" + Utils.realFormat(x, 2) + ") seems quite large!");
    }
  }

}
