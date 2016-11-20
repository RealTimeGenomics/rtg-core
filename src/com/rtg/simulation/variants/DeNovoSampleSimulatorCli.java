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
import com.rtg.util.InvalidParamsException;
import com.rtg.util.PortableRandom;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.intervals.LongRange;
import com.rtg.variant.GenomePriorParams;
import com.rtg.vcf.VcfUtils;

/**
 * Generate a derived genotype that contains additional de novo variants.
 *
 */
public class DeNovoSampleSimulatorCli extends AbstractCli {

  private static final String MODULE_NAME = "denovosim";

  private static final String INPUT_VCF = "input";
  private static final String SAMPLE_FLAG = "sample";
  private static final String ORIGINAL_FLAG = "original";
  private static final String OUTPUT_VCF = "output";
  private static final String OUTPUT_SDF = "output-sdf";
  private static final String EXPECTED_MUTATIONS = "num-mutations";
  private static final String SHOW_MUTATIONS = "show-mutations";
  private static final String REFERENCE_SDF = "reference";
  private static final String SEED = "seed";
  private static final String PRIORS_FLAG = "Xpriors";
  private static final String PLOIDY = "ploidy";

  // We inject one mutation per chromosome, but there is a small probability of an extra mutation.
  private static final int DEFAULT_MUTATIONS_PER_GENOME = 70;


  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return "generate a VCF containing a derived genotype containing de novo variants";
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription("Generates a VCF containing a derived genotype containing de novo variants.");
    mFlags.registerExtendedHelp();
    CommonFlagCategories.setCategories(mFlags);
    mFlags.registerRequired('t', REFERENCE_SDF, File.class, "SDF", "SDF containing input genome").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerRequired('o', OUTPUT_VCF, File.class, "FILE", "output VCF file name").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional(OUTPUT_SDF, File.class, "SDF", "if set, output genome SDF name").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerRequired('i', INPUT_VCF, File.class, "FILE", "input VCF containing genotype of original sample").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerRequired('s', SAMPLE_FLAG, String.class, "STRING", "name for new derived sample").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional('p', PRIORS_FLAG, String.class, "STRING", "selects a properties file specifying the mutation priors. Either a file name or one of [human]", "human").setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerRequired(ORIGINAL_FLAG, String.class, "STRING", "name of the existing sample to use as the original genotype").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional(PLOIDY, ReferencePloidy.class, "string", "ploidy to use", ReferencePloidy.AUTO).setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional(EXPECTED_MUTATIONS, Integer.class, "INT", "expected number of mutations per genome", DEFAULT_MUTATIONS_PER_GENOME).setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional(SEED, Integer.class, "INT", "seed for the random number generator").setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional(SHOW_MUTATIONS, "if set, display information regarding de novo mutation points").setCategory(CommonFlagCategories.UTILITY);
    CommonFlags.initNoGzip(mFlags);

    mFlags.setValidator(new DeNovoSampleSimulatorFlagValidator());
  }

  private static class DeNovoSampleSimulatorFlagValidator implements Validator {

    @Override
    public boolean isValid(final CFlags cflags) {
      if (!cflags.checkNand(OUTPUT_SDF, CommonFlags.NO_GZIP)) {
        return false;
      }
      if (!CommonFlags.validateNotStdout((File) cflags.getValue(OUTPUT_VCF))) {
        return false;
      }
      return !(cflags.isSet(OUTPUT_SDF) && !CommonFlags.validateOutputDirectory(cflags, OUTPUT_SDF));
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
    final String original = (String) flags.getValue(ORIGINAL_FLAG);
    final ReferencePloidy ploidy = (ReferencePloidy) flags.getValue(PLOIDY);
    final GenomePriorParams priors;
    try {
      priors = GenomePriorParams.builder().genomePriors((String) mFlags.getValue(PRIORS_FLAG)).create();
    } catch (final InvalidParamsException e) {
      return 1;
    }
    try (SequencesReader dsr = SequencesReaderFactory.createMemorySequencesReaderCheckEmpty(reference, true, false, LongRange.NONE)) {
      final DeNovoSampleSimulator ss = new DeNovoSampleSimulator(dsr, priors, random, ploidy, (Integer) flags.getValue(EXPECTED_MUTATIONS), flags.isSet(SHOW_MUTATIONS));
      ss.mutateIndividual(popVcf, outputVcf, original, sample);
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
    new DeNovoSampleSimulatorCli().mainExit(args);
  }
}
