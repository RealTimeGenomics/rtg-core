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
package com.rtg.simulation.snpsim;

import static com.rtg.launcher.CommonFlags.NO_GZIP;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.SequencesReader;
import com.rtg.reference.Sex;
import com.rtg.relation.GenomeRelationships;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.PortableRandom;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.Timer;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.GzipAsynchOutputStream;
import com.rtg.util.io.LogStream;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;

/**
 * Command line wrapper for genome mutator (deprecated)
 */
public class GenomeMutatorCli extends LoggedCli {

  private static final String TWIN_DIRECTORY = "diploid";
  private static final String INPUT_DIRECTORY = "input";
  private static final String SNP_FILE = "snp-file";
  private static final String SEED = "seed";
  private static final String PRIORS_FLAG = "Xpriors";

  private static final String MUTATION_RATE = "Xmutation-rate";
  private static final String MUTATION_COUNT = "Xmutation-count";
  private static final String MUTATION_MIN_DISTANCE = "Xmutation-min-distance";
  private static final String SNP_RATE = "Xsnp-rate";
  private static final String MNP_RATE = "Xmnp-rate";
  private static final String INDEL_RATE = "Xindel-rate";
  private static final String VERBOSE = "Xverbose";
  private static final String SIMPLE_MNPS = "Xsimple-mnps";
  private static final String SEX_FLAG = "sex";
  private static final String SAMPLE_FLAG = "sample-name";

  static final double DEFAULT_MUTATION_RATE = 0.001;
  static final int DEFAULT_DISTANCE = 1;

  @TestClass(value = "com.rtg.simulation.snpsim.GenomeMutatorValidatorTest")
  static class GenomeMutatorValidator implements Validator {

    /**
     * Checks if flags are good
     * @param flags the flags
     * @return true if good
     */
    @Override
    public boolean isValid(final CFlags flags) {
      final String inputPath = flags.getFlag(INPUT_DIRECTORY).getValue().toString();
      final File input = new File(inputPath);
      if (!input.exists()) {
        Diagnostic.error(ErrorType.INFO_ERROR, "The specified SDF, \"" + input.getPath() + "\", does not exist.");
        return false;
      }
      if (!input.isDirectory()) {
        Diagnostic.error(ErrorType.INFO_ERROR, "The specified file, \"" + input.getPath() + "\", is not an SDF.");
        return false;
      }

      if (!CommonFlags.validateOutputDirectory(flags)) {
        return false;
      } else if (flags.isSet(TWIN_DIRECTORY) && !CommonFlags.validateOutputDirectory((File) flags.getValue(TWIN_DIRECTORY))) {
        return false;
      }
      final File snpFile = FileUtils.getZippedFileName(!flags.isSet(NO_GZIP), (File) flags.getValue(SNP_FILE));
      if (snpFile.exists()) {
        flags.setParseMessage("SNP mapping file already exists at " + snpFile.toString());
        return false;
      }

      if (flags.isSet(PRIORS_FLAG)) {
        if (flags.isSet(MUTATION_COUNT)) {
          flags.setParseMessage("Cannot set --" + PRIORS_FLAG + " and --" + MUTATION_COUNT);
          return false;
        }
        if (flags.isSet(MUTATION_RATE)) {
          flags.setParseMessage("Cannot set --" + PRIORS_FLAG + " and --" + MUTATION_RATE);
          return false;
        }
        if (flags.isSet(SNP_RATE)) {
          flags.setParseMessage("Cannot set --" + PRIORS_FLAG + " and --" + SNP_RATE);
          return false;
        }
        if (flags.isSet(MNP_RATE)) {
          flags.setParseMessage("Cannot set --" + PRIORS_FLAG + " and --" + MNP_RATE);
          return false;
        }
        if (flags.isSet(INDEL_RATE)) {
          flags.setParseMessage("Cannot set --" + PRIORS_FLAG + " and --" + INDEL_RATE);
          return false;
        }
      }
      // X flags testing
      if (flags.isSet(SEX_FLAG)) {
        if (!flags.isSet(TWIN_DIRECTORY)) {
          flags.setParseMessage("Must specify --" + TWIN_DIRECTORY + " when using --" + SEX_FLAG);
          return false;
        }
        if (Sex.EITHER != flags.getValue(SEX_FLAG) && !CommonFlags.validateSexTemplateReference(flags, SEX_FLAG, null, INPUT_DIRECTORY)) {
          return false;
        }
      }
      if (flags.isSet(MUTATION_COUNT) && flags.isSet(MUTATION_RATE)) {
        flags.setParseMessage("Must specify at most one of --" + MUTATION_COUNT + " and --" + MUTATION_RATE);
        return false;
      }
      if (flags.isSet(MUTATION_RATE)) {
        final Double rate = (Double) flags.getValue(MUTATION_RATE);
        if (rate < 0.0 || rate > 1.0) {
          flags.setParseMessage("Mutation rate must be between 0.0 and 1.0");
          return false;
        }
      }
      if (flags.isSet(MUTATION_COUNT)) {
        final Integer count = (Integer) flags.getValue(MUTATION_COUNT);
        if (count < 0) {
          flags.setParseMessage("Mutation count must be positive");
          return false;
        }
      }
      int numSet = 0;
      double total = 0;
      if (flags.isSet(SNP_RATE)) {
        final double temp = (Double) flags.getValue(SNP_RATE);
        if (temp < 0.0) {
          flags.setParseMessage("SNP rate is under 0");
          return false;
        }
        total += temp;
        numSet++;
      }
      if (flags.isSet(MNP_RATE)) {
        final double temp = (Double) flags.getValue(MNP_RATE);
        if (temp < 0.0) {
          flags.setParseMessage("MNP rate is under 0");
          return false;
        }
        total += temp;
        numSet++;
      }
      if (flags.isSet(INDEL_RATE)) {
        final double temp = (Double) flags.getValue(INDEL_RATE);
        if (temp < 0.0) {
          flags.setParseMessage("INDEL rate is under 0");
          return false;
        }
        total += temp;
        numSet++;
      }
      if (numSet == 3) {
        if (!(Math.abs(total - 1) < 0.00001)) {
          flags.setParseMessage("All (SNP, MNP, INDEL) mutation rates are set but do not total 1.0");
          return false;
        }
      } else if (total > 1.0) {
        flags.setParseMessage("Mutation rates which are set (from SNP, MNP, INDEL) total to greater than 1.0");
        return false;
      }
      return true;
    }
  }

  @Override
  public String moduleName() {
    return "snpsim";
  }

  @Override
  public String description() {
    return "generate a mutated genome by adding SNPs to a reference";
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  /**
   * set up a flags object for this module
   * @param flags the flags to set up
   */
  void initFlags(CFlags flags) {
    flags.registerExtendedHelp();
    flags.setDescription("Creates a simulated mutation of a reference genome. Can generate either haploid or diploid mutations.");
    CommonFlagCategories.setCategories(flags);
    flags.registerRequired('i', INPUT_DIRECTORY, File.class, "SDF", "SDF containing input genome").setCategory(INPUT_OUTPUT);
    flags.registerRequired('o', CommonFlags.OUTPUT_FLAG, File.class, "SDF", "name for genome output SDF").setCategory(INPUT_OUTPUT);
    flags.registerOptional('O', TWIN_DIRECTORY, File.class, "SDF", "name for secondary genome output SDF. Setting this will enable generation of diploid mutations").setCategory(INPUT_OUTPUT);
    flags.registerRequired('s', SNP_FILE, File.class, "FILE", "output file with SNP information").setCategory(INPUT_OUTPUT);
    flags.registerOptional('p', PRIORS_FLAG, String.class, "string", "selects a properties file specifying the priors. Either a file name or one of [human]", "human").setCategory(UTILITY);
    flags.registerOptional(SEED, Integer.class, "INT", "seed for the random number generator").setCategory(UTILITY);

    flags.registerOptional('r', MUTATION_RATE, Double.class, "FLOAT", String.format("fraction of nucleotide positions to be altered (default %1.3f)", DEFAULT_MUTATION_RATE)).setCategory(UTILITY);
    flags.registerOptional('c', MUTATION_COUNT, Integer.class, "INT", "number of mutations to be inserted").setCategory(UTILITY);
    flags.registerOptional('d', MUTATION_MIN_DISTANCE, Integer.class, "INT", "number of nucleotides forced mutation-free between mutations", DEFAULT_DISTANCE).setCategory(UTILITY);
CommonFlags.initNoGzip(flags);
    flags.registerOptional(SNP_RATE, Double.class, "FLOAT", "fraction of mutations that are SNP").setCategory(UTILITY);
    flags.registerOptional(MNP_RATE, Double.class, "FLOAT", "fraction of mutations that are MNP").setCategory(UTILITY);
    flags.registerOptional(INDEL_RATE, Double.class, "FLOAT", "fraction of mutations that are insertions or deletions").setCategory(UTILITY);
    flags.registerOptional(SIMPLE_MNPS, "use old (simple) method to generate MNPs").setCategory(UTILITY);
    flags.registerOptional(VERBOSE, "output into SNP file used for debugging").setCategory(UTILITY);
    flags.registerOptional(SEX_FLAG, Sex.class, "SEX", "sex of individual", Sex.EITHER).setCategory(UTILITY);
    flags.registerOptional(SAMPLE_FLAG, String.class, "STRING", "sample identifier of individual", "SAMPLE").setCategory(UTILITY);

    flags.setValidator(new GenomeMutatorValidator());
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  /**
   * Main method
   * @param args command line arguments
   */
  public static void main(final String[] args) {
    new GenomeMutatorCli().mainExit(args);
  }

  GenomeMutator getMutator(Integer seed, boolean verbose, boolean simpleMnps, int minDistance, VariantParams variantParams, String sample) {
    return new GenomeMutator(seed == null ? new PortableRandom() : new PortableRandom(seed), verbose, simpleMnps, minDistance, variantParams, sample);
  }

  /**
   * Main method with custom streams
   * @param out replacement for standard out
   * @param initLog replacement for standard error
   */
  @Override
  protected int mainExec(OutputStream out, LogStream initLog) throws IOException {
    final CFlags flags = mFlags;
    final GenomeMutator gm;

    final SequenceParams genomeParams = SequenceParams.builder().directory((File) mFlags.getValue(INPUT_DIRECTORY)).mode(SequenceMode.UNIDIRECTIONAL).create();

    final Sex sex = (Sex) flags.getValue(SEX_FLAG);
    final String sample = (String) flags.getValue(SAMPLE_FLAG);
    final GenomeRelationships ped = new GenomeRelationships();
    ped.addGenome(sample, sex);
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    vpb.genome(genomeParams.readerParams());
    vpb.genomeRelationships(ped);
    final VariantParams vp = vpb.create();

    gm = getMutator((Integer) flags.getValue(SEED), flags.isSet(VERBOSE),
        flags.isSet(SIMPLE_MNPS), !flags.isSet(MUTATION_MIN_DISTANCE) ?  2 : (Integer) flags.getValue(MUTATION_MIN_DISTANCE) + 1,
            vp, sample);


    final boolean usePriors = !flags.isSet(SNP_RATE) && !flags.isSet(MNP_RATE) && !flags.isSet(INDEL_RATE) && !flags.isSet(MUTATION_COUNT) && !flags.isSet(MUTATION_RATE);
    if (usePriors) {
      final GenomePriorParams priors = GenomePriorParams.builder()
          .genomePriors((String) mFlags.getValue(PRIORS_FLAG))
          .create();
      gm.setPriors(priors);
    }
    final File outputDirectory = (File) flags.getValue(CommonFlags.OUTPUT_FLAG);
    final File twinDirectory;
    if (flags.isSet(TWIN_DIRECTORY)) {
      twinDirectory = (File) flags.getValue(TWIN_DIRECTORY);
    } else {
      twinDirectory = null;
    }

    final boolean gzip = !flags.isSet(NO_GZIP);
    final File outputFile = FileUtils.getZippedFileName(gzip, (File) flags.getValue(SNP_FILE));

    final int result;
    try (SequencesReader dsr = genomeParams.reader(); OutputStream mappingOutput = FileUtils.createOutputStream(outputFile, gzip, false)) {
      if (usePriors) {
        result = gm.mutatePriors(dsr, sex, outputDirectory, twinDirectory, mappingOutput);
      } else {
        final Double snpRate = flags.isSet(SNP_RATE) ? (Double) flags.getValue(SNP_RATE) : 0;
        final Double mnpRate = flags.isSet(MNP_RATE) ? (Double) flags.getValue(MNP_RATE) : 0;
        final Double indelRate = flags.isSet(INDEL_RATE) ? (Double) flags.getValue(INDEL_RATE) : 0;
        gm.setRates(snpRate, mnpRate, indelRate);

        if (flags.isSet(MUTATION_COUNT)) {
          result = gm.mutateCount(dsr, sex, outputDirectory, twinDirectory, mappingOutput, (Integer) flags.getValue(MUTATION_COUNT));
        } else {
          final double rate = flags.isSet(MUTATION_RATE) ? (Double) flags.getValue(MUTATION_RATE) : DEFAULT_MUTATION_RATE;
          result = gm.mutateRate(dsr, sex, outputDirectory, twinDirectory, mappingOutput, rate);
        }
      }
      gm.summarise(outputDirectory, twinDirectory, out);
    }
    if (gzip && GzipAsynchOutputStream.BGZIP) {
      final Timer indexing = new Timer("SnpIndex");
      indexing.start();
      try {
        new TabixIndexer(outputFile).saveVcfIndex();
      } catch (final UnindexableDataException e) {
        Diagnostic.warning("Cannot produce TABIX index for: " + outputFile + ": " + e.getMessage());
      }
      indexing.stop();
      indexing.log();
    }
    return result;
  }
}
