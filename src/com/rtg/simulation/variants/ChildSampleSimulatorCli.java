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

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.reference.ReferenceGenome.ReferencePloidy;
import com.rtg.reference.Sex;
import com.rtg.util.PortableRandom;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.intervals.LongRange;
import com.rtg.vcf.VcfUtils;

/**
 * Generate the genotypes for a new sample which is the child of two existing samples.
 *
 */
public class ChildSampleSimulatorCli extends AbstractCli {

  private static final String MODULE_NAME = "childsim";

  private static final String INPUT_VCF = "input";
  private static final String SAMPLE_FLAG = "sample";
  private static final String FATHER_FLAG = "father";
  private static final String MOTHER_FLAG = "mother";
  private static final String SEX = "sex";
  private static final String OUTPUT_VCF = "output";
  private static final String OUTPUT_SDF = "output-sdf";
  private static final String EXTRA_CROSSOVERS = "num-crossovers";
  private static final String SHOW_CROSSOVERS = "show-crossovers";
  private static final String REFERENCE_SDF = "reference";
  private static final String SEED = "seed";
  private static final String PLOIDY = "ploidy";

  // We require one crossover, but there is a small probability of an extra crossover.
  private static final double EXTRA_CROSSOVERS_PER_CHROMOSOME = 0.01;


  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return "generate a VCF containing a genotype simulated as a child of two parents";
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription("Generates a VCF containing a genotype simulated as a child of two parents.");
    mFlags.registerExtendedHelp();
    CommonFlagCategories.setCategories(mFlags);
    mFlags.registerRequired('t', REFERENCE_SDF, File.class, "SDF", "SDF containing input genome").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerRequired('o', OUTPUT_VCF, File.class, "FILE", "output VCF file name").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional(OUTPUT_SDF, File.class, "SDF", "if set, output genome SDF name").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerRequired('i', INPUT_VCF, File.class, "FILE", "input VCF containing parent variants").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerRequired('s', SAMPLE_FLAG, String.class, "STRING", "name for new child sample").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerRequired(FATHER_FLAG, String.class, "STRING", "name of the existing sample to use as the father").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerRequired(MOTHER_FLAG, String.class, "STRING", "name of the existing sample to use as the mother").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional(SEX, Sex.class, "SEX", "sex of individual", Sex.EITHER).setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional(PLOIDY, ReferencePloidy.class, "string", "ploidy to use", ReferencePloidy.AUTO).setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional(EXTRA_CROSSOVERS, Double.class, "FLOAT", "likelihood of extra crossovers per chromosome", EXTRA_CROSSOVERS_PER_CHROMOSOME).setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional(SEED, Integer.class, "INT", "seed for the random number generator").setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional(SHOW_CROSSOVERS, "if set, display information regarding haplotype selection and crossover points").setCategory(CommonFlagCategories.UTILITY);
    CommonFlags.initNoGzip(mFlags);

    mFlags.setValidator(new ChildSampleSimulatorValidator());
  }

  private static class ChildSampleSimulatorValidator implements Validator {

    @Override
    public boolean isValid(final CFlags cflags) {
      if (!cflags.checkNand(OUTPUT_SDF, CommonFlags.NO_GZIP)) {
        return false;
      }
      if (!CommonFlags.validateNotStdout((File) cflags.getValue(OUTPUT_VCF))) {
        return false;
      }
      return !(cflags.isSet(OUTPUT_SDF) && !CommonFlags.validateOutputDirectory((File) cflags.getValue(OUTPUT_SDF)));
    }
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

    final File reference = (File) flags.getValue(REFERENCE_SDF);
    final File popVcf = (File) flags.getValue(INPUT_VCF);
    final File outputVcf = VcfUtils.getZippedVcfFileName(!flags.isSet(CommonFlags.NO_GZIP), (File) flags.getValue(OUTPUT_VCF));
    final String sample = (String) flags.getValue(SAMPLE_FLAG);
    final String father = (String) flags.getValue(FATHER_FLAG);
    final String mother = (String) flags.getValue(MOTHER_FLAG);
    final Sex sex = (Sex) flags.getValue(SEX);
    final ReferencePloidy ploidy = (ReferencePloidy) flags.getValue(PLOIDY);
    try (SequencesReader dsr = SequencesReaderFactory.createMemorySequencesReaderCheckEmpty(reference, true, false, LongRange.NONE)) {
      final ChildSampleSimulator ss = new ChildSampleSimulator(dsr, random, ploidy, (Double) flags.getValue(EXTRA_CROSSOVERS), flags.isSet(SHOW_CROSSOVERS));
      ss.mutateIndividual(popVcf, outputVcf, sample, sex, father, mother);
      if (flags.isSet(OUTPUT_SDF)) {
        final SampleReplayer vr = new SampleReplayer(dsr);
        vr.replaySample(outputVcf, (File) flags.getValue(OUTPUT_SDF), sample);
      }
    }
    return 0;
  }

  /**
   * @param args arguments
   */
  public static void main(String[] args) {
    new ChildSampleSimulatorCli().mainExit(args);
  }
}
