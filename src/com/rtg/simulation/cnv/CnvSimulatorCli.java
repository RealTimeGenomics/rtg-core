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
package com.rtg.simulation.cnv;

import static com.rtg.launcher.CommonFlags.NO_GZIP;
import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collection;
import java.util.Locale;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.reader.PrereadType;
import com.rtg.reader.SdfWriter;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.simulation.cnv.CnvSimulator.FixedRegion;
import com.rtg.util.Constants;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.PortableRandom;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogStream;

/**
 * Module wrapper for the standard read type generators
 *
 */
public final class CnvSimulatorCli extends LoggedCli {

  static final String MODULE_NAME = "cnvsim";

  private static final String TWIN_DIRECTORY = "twin";
  private static final String INPUT_DIRECTORY = "input";
  private static final String CNV_FILE = "cnv-file";
  private static final String SEED = "seed";
  private static final String PRIORS_FLAG = "Xpriors";
  private static final String SET_CNV = "Xset-single-cnv";
  private static final String SET_CNV_LENGTH = "Xset-fixed-cnv-length";

  private static final String CNV_COUNT = "cnv-count";
  private static final String CNV_PERCENT = "cnv-percent";

  private static class GenomeMutatorValidator implements Validator {

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
      }
      if (!CommonFlags.validateOutputDirectory(flags, TWIN_DIRECTORY)) {
        return false;
      }

      final File cnvFile = getMappingFile(flags, !flags.isSet(NO_GZIP));
      if (cnvFile.exists()) {
        flags.setParseMessage("CNV file already exists");
        return false;
      }

      if (flags.isSet(SET_CNV)) {
        if (flags.isSet(CNV_COUNT) || flags.isSet(CNV_PERCENT) || flags.isSet(SET_CNV_LENGTH)) {
          flags.setParseMessage("If single CNV regions is set, neither fixed length nor count and percent can be set.");
          return false;
        }

        for (final Object cnvStr : flags.getValues(SET_CNV)) {
          final String[] settings = parseCnvSetting((String) cnvStr);
          if (settings == null) {
            flags.setParseMessage("Bad sytax in CNV setting");
            return false;
          }
        }
      } else if (flags.isSet(SET_CNV_LENGTH)) {
        if (flags.isSet(SET_CNV)) {
          flags.setParseMessage("If fixed CNV regions is set, single CNV cannot be set.");
          return false;
        }

        final String cnvLengthStr = (String) flags.getValue(SET_CNV_LENGTH);
        final String[] settings = parseCnvLengthSetting(cnvLengthStr);
        if (settings == null) {
          flags.setParseMessage("Bad sytax in CNV length setting");
          return false;
        }
      } else if (flags.isSet(CNV_COUNT)) {
        if (flags.isSet(CNV_PERCENT)) {
          flags.setParseMessage("If CNV regions are set with count, percent cannot be set.");
          return false;
        }
        final Integer count = (Integer) flags.getValue(CNV_COUNT);
        if (count < 0) {
          flags.setParseMessage("Number of CNV regions must be positive");
          return false;
        }
      }
      if (flags.isSet(CNV_PERCENT)) {
        final double p = (double) flags.getValue(CNV_PERCENT);
        if (p < 0) {
          flags.setParseMessage("Percentage of CNV region must be positive");
          return false;
        }
        if (p >= 100) {
          flags.setParseMessage("Percentage of CNV region should be smaller 100");
          return false;
        }
      }
      return true;
    }
  }

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return "generate a mutated genome by adding CNVs to a template";
  }

  @Override
  protected void initFlags() {
    mFlags.registerExtendedHelp();
    mFlags.setDescription("Creates copy number variation between two simulated genomes, one the baseline and the other the sample to be tested.");
    CommonFlagCategories.setCategories(mFlags);
    mFlags.registerRequired('i', INPUT_DIRECTORY, File.class, "SDF", "SDF containing input genome").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('o', OUTPUT_FLAG, File.class, "SDF", "name for genome output SDF").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('O', TWIN_DIRECTORY, File.class, "SDF", "name for secondary output SDF").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('s', CNV_FILE, File.class, "FILE", "output file with CNV information").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional('p', PRIORS_FLAG, String.class, "string", "selects a properties file specifying the CNV priors. Either a file name or one of ", "cnv-default").setCategory(UTILITY);
    mFlags.registerOptional('n', CNV_COUNT, Integer.class, "INT", "number of regions generated with CNV effects on").setCategory(UTILITY);
    mFlags.registerOptional('N', CNV_PERCENT, Integer.class, "INT", "approximate minimum percent of template region generated with CNV effects on", 10).setCategory(UTILITY);
    CommonFlags.initNoGzip(mFlags);
    Flag flag = mFlags.registerOptional('C', SET_CNV, String.class, "string", "set one CNV; sequence-name:start:end:copy-number");
    flag.setMaxCount(Integer.MAX_VALUE);
    flag.setCategory(UTILITY);
    flag = mFlags.registerOptional('L', SET_CNV_LENGTH, String.class, "string", "set fixed CNV length with max bounds for derivation; length:max-derivation");
    flag.setMaxCount(Integer.MAX_VALUE);
    flag.setCategory(UTILITY);
    mFlags.registerOptional(SEED, Integer.class, "INT", "seed for the random number generator").setCategory(UTILITY);
    mFlags.setValidator(new GenomeMutatorValidator());
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(OUTPUT_FLAG);
  }

  /**
   * Main method
   * @param args command line arguments
   */
  public static void main(final String[] args) {
    new CnvSimulatorCli().mainExit(args);
  }

  /**
   * Main method with custom streams
   * @param out replacement for standard out
   * @param initLog replacement for standard error
   */
  @Override
  protected int mainExec(OutputStream out, LogStream initLog) throws IOException {
    final CFlags flags = mFlags;
    final CnvSimulator simulator;

    final File input = (File) flags.getValue(INPUT_DIRECTORY);
    final File outputDirectory = (File) flags.getValue(OUTPUT_FLAG);
    final File twinDirectory = (File) flags.getValue(TWIN_DIRECTORY);
    final boolean gzip = !flags.isSet(NO_GZIP);

    try {
      try (OutputStream mappingOutput = FileUtils.createOutputStream(getMappingFile(flags, gzip), gzip, false); SequencesReader dsr = SequencesReaderFactory.createDefaultSequencesReaderCheckEmpty(input)) {
        final PrereadType pType = dsr.getPrereadType();

        try (SdfWriter output = new SdfWriter(outputDirectory, Constants.MAX_FILE_SIZE, pType, false, true, !flags.isSet(NO_GZIP), dsr.type());
             SdfWriter twin = new SdfWriter(twinDirectory, Constants.MAX_FILE_SIZE, pType, false, true, !flags.isSet(NO_GZIP), dsr.type())) {

          final PortableRandom random;
          if (mFlags.isSet(SEED)) {
            random = new PortableRandom((Integer) mFlags.getValue(SEED));
          } else {
            random = new PortableRandom();
          }
          final CnvPriorParams priors = CnvPriorParams.builder()
            .cnvpriors((String) mFlags.getValue(PRIORS_FLAG))
            .create();
          if (mFlags.isSet(SET_CNV_LENGTH)) {
            simulator = new CnvSimulator(dsr, output, twin, mappingOutput, random, priors,
              (double) mFlags.getValue(CNV_PERCENT),
              (mFlags.isSet(CNV_COUNT)) ? (int) mFlags.getValue(CNV_COUNT) : Integer.MAX_VALUE);
            simulator.generate(parseCnvLengthSetting((String) mFlags.getValue(SET_CNV_LENGTH)));
          } else if (mFlags.isSet(SET_CNV)) {
            simulator = new CnvSimulator(dsr, output, twin, mappingOutput, random, priors,
              -1.0,
              mFlags.getValues(SET_CNV).size());
            final FixedRegion[] fixedRegions = getRegionArray(mFlags.getValues(SET_CNV));
            simulator.generate(fixedRegions);
          } else {
            simulator = new CnvSimulator(dsr, output, twin, mappingOutput, random, priors,
              (double) mFlags.getValue(CNV_PERCENT),
              (mFlags.isSet(CNV_COUNT)) ? (int) mFlags.getValue(CNV_COUNT) : Integer.MAX_VALUE);
            simulator.generate();
          }
          simulator.outputHistograms(out);
          return 0;
        }
      }
    } catch (final FileNotFoundException e) {
      throw new NoTalkbackSlimException(ErrorType.SDF_INDEX_NOT_VALID, input.toString());
    }
  }

  private static File getMappingFile(final CFlags flags, final boolean gzip) {
    final File mappingFile = (File) flags.getValue(CNV_FILE);
    final boolean zipSuff = gzip && !mappingFile.getAbsolutePath().toLowerCase(Locale.getDefault()).endsWith(FileUtils.GZ_SUFFIX);
    return zipSuff ? new File(mappingFile.getAbsolutePath() + FileUtils.GZ_SUFFIX) : mappingFile;
  }

  private FixedRegion[] getRegionArray(Collection<Object> arrayList) {
    final FixedRegion[] regions = new FixedRegion[arrayList.size()];
    int i = 0;
    for (final Object reg : arrayList) {
      final String[] settings = parseCnvSetting((String) reg);
      if (settings == null) {
        throw new InvalidParamsException(ErrorType.INVALID_PARAMETER, "Syntax error in values");
      }
      regions[i] = new FixedRegion(settings[0], Integer.parseInt(settings[1]), Integer.parseInt(settings[2]),
          Integer.parseInt(settings[3]));
      ++i;
    }
    return regions;
  }

  private static String[] parseCnvSetting(String str) {
    final String[] parts = str.split(":");
    if (parts.length != 4) {
      return null;
    }
    try {
      for (int i = 1; i < parts.length; ++i) {
        final int num = Integer.parseInt(parts[i]);
        if (num < 0) {
          return null;
        }
      }
    } catch (final NumberFormatException e) {
      return null;
    }
    return parts;
  }

  private static String[] parseCnvLengthSetting(String str) {
    final String[] parts = str.split(":");
    if (parts.length != 2) {
      return null;
    }
    try {
      for (final String part : parts) {
        final int num = Integer.parseInt(part);
        if (num < 0) {
          return null;
        }
      }
    } catch (final NumberFormatException e) {
      return null;
    }
    return parts;
  }
}
