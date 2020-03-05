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

import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collection;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.OutputParams;
import com.rtg.launcher.ParamsCli;
import com.rtg.launcher.SequenceParams;
import com.rtg.sam.SamFilterOptions;
import com.rtg.util.IORunnable;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Perform a species analysis of a metagenomics sample.
 */
public class SpeciesCli extends ParamsCli<SpeciesParams> {

  static final char COMMENT_CHAR = '#';
  /** Command line flag for the species genome database */
  public static final String TEMPLATE_FLAG = "genomes";
  private static final String PRINT_ALL_FLAG = "print-all";
  private static final String MIN_CONFIDENCE_VALUE_FLAG = "min-confidence";
  private static final String RELABEL_SPECIES_FLAG = "relabel-species-file";
  private static final String ITERATIONS_FLAG = "Xiterations";
  private static final String STD_DEV_GRAPH_FLAG = "Xstddev-graph";
  private static final String VERBOSE_FLAG = "Xverbose";
  private static final String NAMESPACE_FLAG = "Xsimple-names";

  private static class SpeciesFlagsValidator implements Validator {
    @Override
    public boolean isValid(final CFlags flags) {
      if (flags.isSet(RELABEL_SPECIES_FLAG)) {
        final File referenceMap = (File) flags.getValue(RELABEL_SPECIES_FLAG);
        if (!referenceMap.exists()) {
          flags.setParseMessage("Given input file \"" + referenceMap.getPath() + "\" for --" + RELABEL_SPECIES_FLAG + " does not exist");
          return false;
        } else if (referenceMap.isDirectory()) {
          flags.setParseMessage("Given input file \"" + referenceMap.getPath() + "\" for --" + RELABEL_SPECIES_FLAG + " is a directory, must be a file");
          return false;
        }
      }
      return CommonFlags.validateThreads(flags)
        && CommonFlags.validateSDF(flags, TEMPLATE_FLAG)
        && CommonFlags.validateOutputDirectory(flags)
        && flags.checkInRange(MIN_CONFIDENCE_VALUE_FLAG, 0.0, Double.MAX_VALUE)
        && CommonFlags.checkFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, Integer.MAX_VALUE)
        && SamFilterOptions.validateFilterFlags(flags, false);
    }
  }

  @Override
  public String moduleName() {
    return "species";
  }

  @Override
  public String description() {
    return "estimate species frequency in metagenomic samples";
  }

  @Override
  protected SpeciesParams makeParams() throws IOException {
    final Collection<File> inputFiles = CommonFlags.getFileList(mFlags, CommonFlags.INPUT_LIST_FLAG, null, false);
    Diagnostic.userLog("Input SAM files: " + inputFiles);
    final File output = (File) mFlags.getValue(OUTPUT_FLAG);
    final OutputParams outParams = new OutputParams(output, false);
    final File genomes = (File) mFlags.getValue(TEMPLATE_FLAG);
    final SequenceParams genomesParams = SequenceParams.builder().directory(genomes).create();
    final int minIter = (Integer) mFlags.getValue(ITERATIONS_FLAG);
    final File referenceMap = mFlags.isSet(RELABEL_SPECIES_FLAG) ? (File) mFlags.getValue(RELABEL_SPECIES_FLAG) : null;
    final double minConfidence = (Double) mFlags.getValue(MIN_CONFIDENCE_VALUE_FLAG);
    return SpeciesParams.builder()
      .outputParams(outParams)
      .mapped(inputFiles)
      .filterParams(SamFilterOptions.makeFilterParamsBuilder(mFlags).create())
      .genome(genomesParams.readerParams())
      .minIter(minIter)
      .verbose(mFlags.isSet(VERBOSE_FLAG))
      .referenceMap(referenceMap)
      .printAll(mFlags.isSet(PRINT_ALL_FLAG))
      .execThreads(CommonFlags.parseThreads((Integer) mFlags.getValue(CommonFlags.THREADS_FLAG)))
      .minConfidence(minConfidence)
      .identifierCreator(mFlags.isSet(NAMESPACE_FLAG) ? new DefaultIdentifierCreator() : new HashingIdentifierCreator())
      .create();
  }

  @Override
  protected IORunnable task(final SpeciesParams params, final OutputStream out) {
    return new SpeciesTask(params, out, mUsageMetric);
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
    CommonFlagCategories.setCategories(flags);
    flags.registerRequired('o', OUTPUT_FLAG, File.class, CommonFlags.DIR, "directory for output").setCategory(INPUT_OUTPUT);
    flags.registerRequired('t', TEMPLATE_FLAG, File.class, CommonFlags.SDF, "SDF containing the genomes").setCategory(INPUT_OUTPUT);
    flags.registerOptional('r', RELABEL_SPECIES_FLAG, File.class, CommonFlags.FILE, "file containing list of species name to reference name mappings (1 mapping per line format: [reference short name][tab][species])").setCategory(INPUT_OUTPUT);
    flags.registerOptional(ITERATIONS_FLAG, Integer.class, CommonFlags.INT, "minimum number of iterations multiplied by the block size", 1).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(PRINT_ALL_FLAG, "print non present species in the output file").setCategory(REPORTING);
    flags.registerOptional('c', MIN_CONFIDENCE_VALUE_FLAG, Double.class, CommonFlags.FLOAT, "species below this confidence value will not be reported", 10.0).setCategory(REPORTING);
    flags.registerOptional(NAMESPACE_FLAG, "use only the read name to identify sequences").setCategory(UTILITY);

    CommonFlags.initThreadsFlag(flags);
    final Flag<File> listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, CommonFlags.FILE, "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);
    flags.registerOptional(VERBOSE_FLAG, "turn on output of convergence information").setCategory(UTILITY);
    flags.registerOptional(STD_DEV_GRAPH_FLAG, "output graph across 1 standard deviation").setCategory(INPUT_OUTPUT);
    //You should never filter in IH (it is expected that there will be hits on multiple bacterial species)
    //SamFilterOptions.registerMaxHitsFlag(flags, 'c');
    SamFilterOptions.registerMaxASMatedFlag(flags, 'm');
    SamFilterOptions.registerMaxASUnmatedFlag(flags, 'u');
    SamFilterOptions.registerExcludeMatedFlag(flags);
    SamFilterOptions.registerExcludeUnmatedFlag(flags);
    final Flag<File> inFlag = flags.registerRequired(File.class, CommonFlags.FILE, "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    flags.addRequiredSet(inFlag);
    flags.addRequiredSet(listFlag);
  }

}

