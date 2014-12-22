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
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.PortableRandom;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;

/**
 * Command line wrapper for fixed step variant generator
 */
public class FixedStepPopulationVariantGeneratorCli extends AbstractCli {

  private static final String MODULE_NAME = "fixedstepsnpsim";
  private static final String OUTPUT_VCF = "output";
  private static final String REFERENCE_SDF = "reference";
  private static final String DISTANCE = "distance";
  private static final String SEED = "seed";
  private static final String SNP_SPECIFICATION = "spec";
  private static final String FREQUENCY = "allele-frequency";

  /**
   * @return current name of the module
   */
  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  /**
   * Set up a flags object for this module
   * @param flags the flags to set up
   */
  public void initFlags(CFlags flags) {
    flags.registerExtendedHelp();
    CommonFlagCategories.setCategories(flags);
    flags.registerRequired('i', REFERENCE_SDF, File.class, "SDF", "SDF containing input genome").setCategory(INPUT_OUTPUT);
    flags.registerRequired('o', OUTPUT_VCF, File.class, "FILE", "name for population output VCF").setCategory(INPUT_OUTPUT);
    flags.registerRequired(SNP_SPECIFICATION, String.class, "string", "generated mutation format").setCategory(INPUT_OUTPUT);
    flags.registerOptional(SEED, Integer.class, "INT", "seed for the random number generator").setCategory(UTILITY);
    flags.registerRequired('d', DISTANCE, Integer.class, "INT", "distance between mutations").setCategory(INPUT_OUTPUT);
    flags.registerOptional('a', FREQUENCY, Double.class, "FLOAT", "allele frequency", 0.5).setCategory(UTILITY);
    flags.setValidator(new Validator() {
      @Override
      public boolean isValid(CFlags flags) {
        final Double af = (Double) flags.getValue(FREQUENCY);
        return af != null && af >= 0 && af <= 1 && !af.isNaN();
      }
    });
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

    final int distance = (Integer) flags.getValue(DISTANCE);
    final File input = (File) flags.getValue(REFERENCE_SDF);
    final Mutator mutator = new Mutator((String) flags.getValue(SNP_SPECIFICATION));
    final File outputVcf = FileUtils.getZippedFileName(true, (File) flags.getValue(OUTPUT_VCF));
    final double af = (Double) flags.getValue(FREQUENCY);
    try (SequencesReader dsr = SequencesReaderFactory.createMemorySequencesReader(input, true, LongRange.NONE)) {
      final FixedStepPopulationVariantGenerator fs = new FixedStepPopulationVariantGenerator(dsr, distance, mutator, random, af);
      PopulationVariantGenerator.writeAsVcf(outputVcf, null, fs.generatePopulation(), dsr);
      return 0;
    }
  }

  /**
   * Main method
   * @param args command line arguments
   */
  public static void main(final String[] args) {
    new FixedStepPopulationVariantGeneratorCli().mainExit(args);
  }

}
