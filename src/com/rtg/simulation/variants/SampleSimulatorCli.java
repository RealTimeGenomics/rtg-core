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
 * Randomly generate the genotypes for a new sample, given a set of population variants
 * with allele frequencies.
 *
 */
public class SampleSimulatorCli extends AbstractCli {

  private static final String MODULE_NAME = "samplesim";

  private static final String POPULATION_VCF = "input";
  private static final String SAMPLE_NAME = "sample";
  private static final String SEX = "sex";
  private static final String OUTPUT_VCF = "output";
  private static final String OUTPUT_SDF = "output-sdf";
  private static final String REFERENCE_SDF = "reference";
  private static final String SEED = "seed";
  private static final String PLOIDY = "ploidy";


  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return "generate a VCF containing a genotype simulated from a population";
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription("Generates a VCF containing a genotype simulated from a population.");
    mFlags.registerExtendedHelp();
    CommonFlagCategories.setCategories(mFlags);
    mFlags.registerRequired('t', REFERENCE_SDF, File.class, "SDF", "SDF containing input genome").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerRequired('o', OUTPUT_VCF, File.class, "FILE", "output VCF file name").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional(OUTPUT_SDF, File.class, "SDF", "if set, output genome SDF name").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerRequired('i', POPULATION_VCF, File.class, "FILE", "input VCF containing population variants").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerRequired('s', SAMPLE_NAME, String.class, "STRING", "name for sample").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional(SEX, Sex.class, "SEX", "sex of individual", Sex.EITHER).setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional(PLOIDY, ReferencePloidy.class, "string", "ploidy to use", ReferencePloidy.AUTO).setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional(SEED, Integer.class, "INT", "seed for the random number generator").setCategory(CommonFlagCategories.UTILITY);
    CommonFlags.initNoGzip(mFlags);

    mFlags.setValidator(new SampleSimulatorFlagValidator());
  }

  private static class SampleSimulatorFlagValidator implements Validator {

    @Override
    public boolean isValid(final CFlags cflags) {
      if (cflags.isSet(OUTPUT_SDF) && !CommonFlags.validateOutputDirectory(cflags, OUTPUT_SDF)) {
        return false;
      }
      if (!cflags.checkNand(OUTPUT_SDF, CommonFlags.NO_GZIP)) {
        return false;
      }
      return CommonFlags.validateNotStdout((File) cflags.getValue(OUTPUT_VCF));
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
    final File popVcf = (File) flags.getValue(POPULATION_VCF);
    final File outputVcf = VcfUtils.getZippedVcfFileName(!flags.isSet(CommonFlags.NO_GZIP), (File) flags.getValue(OUTPUT_VCF));
    final String sample = (String) flags.getValue(SAMPLE_NAME);
    final Sex sex = (Sex) flags.getValue(SEX);
    final ReferencePloidy ploidy = (ReferencePloidy) flags.getValue(PLOIDY);
    try (SequencesReader dsr = SequencesReaderFactory.createMemorySequencesReaderCheckEmpty(reference, true, false, LongRange.NONE)) {
      final SampleSimulator ss = new SampleSimulator(dsr, random, ploidy);
      ss.mutateIndividual(popVcf, outputVcf, sample, sex);
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
    new SampleSimulatorCli().mainExit(args);
  }
}
