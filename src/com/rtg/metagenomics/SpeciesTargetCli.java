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
package com.rtg.metagenomics;

import static com.rtg.launcher.BuildCommon.RESOURCE;
import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collection;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.OutputParams;
import com.rtg.launcher.ParamsCli;
import com.rtg.launcher.SequenceParams;
import com.rtg.sam.SamFilterOptions;
import com.rtg.util.Constants;
import com.rtg.util.IORunnable;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Perform a species target analysis of a metagenomics sample.
 */
public class SpeciesTargetCli extends ParamsCli<SpeciesParams> {

  private static final String TARGET_FLAG = "target";
  /** Command line flag for the species genome database */
  public static final String TEMPLATE_FLAG = "genomes";

  private static class SpeciesFlagsValidator implements Validator {
    @Override
    public boolean isValid(final CFlags flags) {
      if (!CommonFlags.validateOutputDirectory(flags)) {
        return false;
      }
      if (!CommonFlags.validateSDF(flags, TEMPLATE_FLAG)) {
        return false;
      }

      if (!CommonFlags.validateThreads(flags)) {
        return false;
      }

      if (!CommonFlags.checkFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, Constants.MAX_OPEN_FILES)) {
        return false;
      }

      if (!SamFilterOptions.validateFilterFlags(flags, false)) {
        return false;
      }

      return true;
    }
  }

  @Override
  protected SpeciesParams makeParams() throws InvalidParamsException, IOException {
    final Collection<File> inputFiles = CommonFlags.getFileList(mFlags, CommonFlags.INPUT_LIST_FLAG, null, false);
    Diagnostic.userLog("Input SAM files: " + inputFiles);
    final File output = (File) mFlags.getValue(OUTPUT_FLAG);
    final OutputParams outParams = new OutputParams(output, false, false);
    final File genomes = (File) mFlags.getValue(TEMPLATE_FLAG);
    final SequenceParams genomesParams = SequenceParams.builder().directory(genomes).create();
    return SpeciesParams.builder()
        .outputParams(outParams)
        .mapped(inputFiles)
        .filterParams(SamFilterOptions.makeFilterParamsBuilder(mFlags).excludeUnmapped(true).excludeUnplaced(true).create())
        .genome(genomesParams)
        .execThreads(CommonFlags.parseThreads((Integer) mFlags.getValue(CommonFlags.THREADS_FLAG)))
        .create();
  }

  @Override
  protected IORunnable task(final SpeciesParams params, final OutputStream out) {
    return new SpeciesTargetTask(params, out, mUsageMetric, (Integer) mFlags.getValue(TARGET_FLAG));
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(OUTPUT_FLAG);
  }

  @Override
  protected void initFlags() {
    initFlagsLocal(mFlags);
  }

  /**
   * Initialise a set of flags for the species module
   * @param flags the flags object to initialise
   */
  public static void initFlagsLocal(CFlags flags) {
    flags.setValidator(new SpeciesFlagsValidator());
    flags.setDescription("Calculates a species distribution from a metagenomic sample.");
    flags.registerExtendedHelp();
    CommonFlagCategories.setCategories(flags);
    flags.registerRequired('o', OUTPUT_FLAG, File.class, "DIR", RESOURCE.getString("OUTPUT_DESC")).setCategory(INPUT_OUTPUT);
    flags.registerRequired('t', TEMPLATE_FLAG, File.class, "SDF", "SDF containing the genomes").setCategory(INPUT_OUTPUT);
    flags.registerRequired(TARGET_FLAG, Integer.class, "INT", "identifier for target species").setCategory(INPUT_OUTPUT);

    CommonFlags.initThreadsFlag(flags);
    final Flag listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, "FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);
    SamFilterOptions.registerMaxASMatedFlag(flags, 'm');
    SamFilterOptions.registerMaxASUnmatedFlag(flags, 'u');
    SamFilterOptions.registerExcludeMatedFlag(flags);
    SamFilterOptions.registerExcludeUnmatedFlag(flags);
    final Flag inFlag = flags.registerRequired(File.class, "FILE", "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    flags.addRequiredSet(inFlag);
    flags.addRequiredSet(listFlag);
  }

  @Override
  public String moduleName() {
    return "species-target";
  }

  /**
   * Main program for species metagenomics. Use -h to get help.
   * @param args command line arguments.
   */
  public static void main(final String[] args) {
    new SpeciesTargetCli().mainExit(args);
  }

}

