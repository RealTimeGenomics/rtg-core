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
package com.rtg.simulation.variants;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.PortableRandom;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;
import com.rtg.variant.GenomePriorParams;

/**
 * Command line wrapper for prior-based population variant creation
 */
public class PriorPopulationVariantGeneratorCli extends AbstractCli {

  private static final String MODULE_NAME = "popsim";
  private static final String OUTPUT_VCF = "output";
  private static final String REFERENCE_SDF = "reference";
  private static final String SEED = "seed";
  private static final String PRIORS_FLAG = "Xpriors";
  private static final String BIAS_FLAG = "Xbias";
  private static final String RATE_FLAG = "Xrate";

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return "generate a VCF containing simulated population variants";
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  /**
   * set up a flags object for this module
   *
   * @param flags the flags to set up
   */
  public void initFlags(CFlags flags) {
    flags.setDescription("Generates a VCF containing simulated population variants.");
    flags.registerExtendedHelp();
    CommonFlagCategories.setCategories(flags);
    flags.registerRequired('t', REFERENCE_SDF, File.class, "SDF", "SDF containing input genome").setCategory(INPUT_OUTPUT);
    flags.registerRequired('o', OUTPUT_VCF, File.class, "FILE", "output VCF file name").setCategory(INPUT_OUTPUT);
    flags.registerOptional('p', PRIORS_FLAG, String.class, "STRING", "selects a properties file specifying the priors. Either a file name or one of [human]", "human").setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(BIAS_FLAG, Double.class, "FLOAT", "bias frequency of variants towards alt alleles.", 0.0).setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(RATE_FLAG, Double.class, "FLOAT", "per base rate of variant generation (overrides that loaded from priors).").setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(SEED, Integer.class, "INT", "seed for the random number generator").setCategory(UTILITY);
    CommonFlags.initNoGzip(flags);
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final CFlags flags = mFlags;
    final PortableRandom random;
    if (flags.isSet(SEED)) {
      random = new PortableRandom((Integer) flags.getValue(SEED));
    } else {
      random = new PortableRandom();
    }

    final PopulationMutatorPriors priors;
    try {
      priors = new PopulationMutatorPriors(GenomePriorParams.builder().genomePriors((String) mFlags.getValue(PRIORS_FLAG)).create());
    } catch (final InvalidParamsException e) {
      return 1;
    }
    final File reference = (File) flags.getValue(REFERENCE_SDF);
    final File outputVcf = FileUtils.getZippedFileName(!flags.isSet(CommonFlags.NO_GZIP), (File) flags.getValue(OUTPUT_VCF));
    try (SequencesReader dsr = SequencesReaderFactory.createMemorySequencesReaderCheckEmpty(reference, true, false, LongRange.NONE)) {
      final int targetVariants;
      if (flags.isSet(RATE_FLAG)) {
        targetVariants = (int) (dsr.totalLength() * (Double) flags.getValue(RATE_FLAG));
      } else {
        targetVariants = (int) (dsr.totalLength() * priors.rate());
      }
      final PriorPopulationVariantGenerator fs = new PriorPopulationVariantGenerator(dsr, priors, random, (Double) flags.getValue(BIAS_FLAG), targetVariants);
      PopulationVariantGenerator.writeAsVcf(outputVcf, null, fs.generatePopulation(), dsr);
    }
    return 0;
  }

  /**
   * Main method
   *
   * @param args command line arguments
   */
  public static void main(final String[] args) {
    new PriorPopulationVariantGeneratorCli().mainExit(args);
  }

}
