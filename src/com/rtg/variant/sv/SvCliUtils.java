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

package com.rtg.variant.sv;

import static com.rtg.launcher.BuildCommon.RESOURCE;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import com.rtg.launcher.BuildCommon;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.OutputParams;
import com.rtg.launcher.SequenceParams.SequenceParamsBuilder;
import com.rtg.sam.SamFilterOptions;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
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
      if (!CommonFlags.validateOutputDirectory(flags)) {
        return false;
      }
      if (!CommonFlags.validateSDF(flags, CommonFlags.TEMPLATE_FLAG)) {
        return false;
      }
      if (!CommonFlags.checkFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, Integer.MAX_VALUE)) {
        return false;
      }
      if (!CommonFlags.checkFileList(flags, RG_STATS_LIST_FLAG, RG_STATS_FILE, Integer.MAX_VALUE)) {
        return false;
      }
      if (!CommonFlags.validateThreads(flags)) {
        return false;
      }
      if (!SamFilterOptions.validateFilterFlags(flags, false)) {
        return false;
      }
      return true;
    }
  }

  /**
   * Initialise the flags common to all structural variant commands
   * @param flags the flags to initialise
   */
  public static void initCommonFlags(CFlags flags) {
    flags.registerExtendedHelp();
    CommonFlagCategories.setCategories(flags);
    flags.registerRequired('t', CommonFlags.TEMPLATE_FLAG, File.class, CommonFlags.SDF, "SDF of the reference genome the reads have been mapped against").setCategory(INPUT_OUTPUT);
    final Flag inFlag = flags.registerRequired(File.class, "FILE", "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    final Flag listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, "FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);
    flags.registerRequired('o', CommonFlags.OUTPUT_FLAG, File.class, "DIR", RESOURCE.getString("OUTPUT_DESC")).setCategory(INPUT_OUTPUT);
    CommonFlags.initNoGzip(flags);
    final Flag rgstats = flags.registerOptional('r', RG_STATS_FILE, File.class, "FILE", "text file containing read group stats").setCategory(INPUT_OUTPUT);
    rgstats.setMinCount(0);
    rgstats.setMaxCount(Integer.MAX_VALUE);
    final Flag rgListFlag = flags.registerOptional('R', RG_STATS_LIST_FLAG, File.class, "FILE", "file containing list of read group stats files (1 per line)").setCategory(INPUT_OUTPUT);
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
    flags.registerOptional(RELABEL_FLAG, File.class, "FILE", "file containing read group relabel mappings (1 per line)").setCategory(INPUT_OUTPUT);
  }

  /**
   * Populate a builder with the common flag settings
   * @param builder the builder to populate
   * @param genomeBuilder <code>SequenceParamsBuilder</code> to add the directory to and use
   * @param flags the flags to use to populate the builder
   * @throws InvalidParamsException if there are errors in the values of the command line flags
   * @throws IOException If an I/O error occurs
   */
  public static void populateCommonParams(SvParamsBuilder<?> builder, SequenceParamsBuilder genomeBuilder, CFlags flags)  throws InvalidParamsException, IOException {
    builder.name(flags.getName())
           .genome(genomeBuilder.directory((File) flags.getValue(CommonFlags.TEMPLATE_FLAG)).create())
           .outputParams(new OutputParams((File) flags.getValue(CommonFlags.OUTPUT_FLAG), flags.isSet(BuildCommon.PROGRESS_FLAG), !flags.isSet(CommonFlags.NO_GZIP)))
           .ioThreads(CommonFlags.parseIOThreads((Integer) flags.getValue(CommonFlags.THREADS_FLAG)))
           .mapped(CommonFlags.getFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, false))
           .filterParams(SamFilterOptions.makeFilterParamsBuilder(flags).excludeUnmapped(true).excludeUnplaced(true).create());

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
    final Map<String, ReadGroupStats> readGroupStats = ReadGroupStats.loadReadGroupStats(relabelMap, statsFiles.toArray(new File[statsFiles.size()]));
    for (final ReadGroupStats e : readGroupStats.values()) {
      if (!e.isValid()) {
        throw new NoTalkbackSlimException("Invalid read group statistics for read group: " + e.id());
      //the below checks are what was previously being asserted upon output of rg stats file.
      } else if (e.gapMean() <= 20.0) {
        throw new NoTalkbackSlimException("Gap mean in statistics for read group: " + e.id() + " must be greater than 20.0");
      } else if (e.meanLength() <= 10.0) {
        throw new NoTalkbackSlimException("Mean length in statistics for read group: " + e.id() + " must be greater than 10.0");
      } else if (e.fragmentMean() <= 20.0) {
        throw new NoTalkbackSlimException("Fragment mean in statistics for read group: " + e.id() + " must be greater than 20.0");
      } else if (e.fragmentStdDev() == 0.0 || e.gapStdDev() == 0.0 || e.properRate() == 0.0 || e.properRandomRate() == 0.0 || e.discordantRate() == 0.0 || e.unmatedRate() == 0.0) {
        throw new NoTalkbackSlimException("Statistics for read group: " + e.id() + " must be greater than 0.0");
      }
    }
    builder.readGroupStatistics(readGroupStats);
  }

}
