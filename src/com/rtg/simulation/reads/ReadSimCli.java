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
package com.rtg.simulation.reads;

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.bed.BedUtils;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.mode.SequenceType;
import com.rtg.reader.PrereadNamesInterface;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.sam.SamCommandHelper;
import com.rtg.simulation.genome.SequenceDistribution;
import com.rtg.taxonomy.TaxonomyUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.MathUtils;
import com.rtg.util.PortableRandom;
import com.rtg.util.Utils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.WarningType;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogStream;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.MachineErrorParams;
import com.rtg.variant.MachineErrorParamsBuilder;

/**
 * Module wrapper for the standard read type generators
 */
public class ReadSimCli extends LoggedCli {

  static final String MODULE_NAME = "readsim";

  // Common to all machines
  static final String INPUT = "input";
  static final String TWIN_INPUT = "Xdiploid-input";
  static final String OUTPUT = "output";
  static final String SEED = "seed";
  static final String COMMENT = "comment";
  static final String MACHINE_TYPE = "machine";
  static final String MACHINE_ERROR_PRIORS = "Xmachine-errors";
  static final String COVERAGE = "coverage";
  static final String READS = "num-reads";
  static final String DISTRIBUTION = "distribution";
  static final String TAXONOMY_DISTRIBUTION = "taxonomy-distribution";
  static final String ABUNDANCE = "abundance";
  static final String DNA_FRACTION = "dna-fraction";
  static final String N_RATE = "n-rate";
  static final String QUAL_RANGE = "qual-range";

  // Options for library generation
  static final String MAX_FRAGMENT = "max-fragment-size";   // Fragment size from which reads are taken
  static final String MIN_FRAGMENT = "min-fragment-size";
  static final String ALLOW_UNKNOWNS = "allow-unknowns";
  //static final String TRIM_FRAGMENT_BELL = "Xtrim-fragment-bell"; // Not currently implemented

  // Single-end Illumina machine
  static final String READLENGTH = "read-length";

  // Paired-end Illumina machine
  static final String LEFT_READLENGTH = "left-read-length";
  static final String RIGHT_READLENGTH = "right-read-length";

  // 454 machine
  static final String MIN_TOTAL_454_LENGTH = "454-min-total-size";
  static final String MAX_TOTAL_454_LENGTH = "454-max-total-size";

  // iontorrent machine
  static final String MIN_TOTAL_IONTORRENT_LENGTH = "ion-min-total-size";
  static final String MAX_TOTAL_IONTORRENT_LENGTH = "ion-max-total-size";

  static final String CAT_FRAGMENTS = "Fragment Generation";
  private static final String CAT_ILLUMINA_SE = "Illumina SE";
  private static final String CAT_ILLUMINA_PE = "Illumina PE";
  private static final String CAT_454_PE = "454 SE/PE";
  private static final String CAT_ION_SE = "IonTorrent SE";

  private static final String NO_NAMES = "no-names";
  private static final String NO_QUAL = "no-qualities";

  private static final String MNP_EVENT_RATE = "Xmnp-event-rate";
  private static final String INS_EVENT_RATE = "Xinsert-event-rate";
  private static final String DEL_EVENT_RATE = "Xdelete-event-rate";
  private static final String BED_FILE = "Xbed-file";

  private static final String PCR_DUP_RATE = "Xpcr-duplicate-rate";
  private static final String CHIMERA_RATE = "Xchimera-rate";


  private static boolean okProbability(final Object o) {
    if (!(o instanceof Double)) {
      return false;
    }
    final double d = (Double) o;
    return d >= 0 && d <= 1;
  }

  @TestClass("com.rtg.simulation.reads.ReadSimValidatorTest")
  protected static class ReadSimValidator implements Validator {

    @Override
    public boolean isValid(final CFlags cflags) {

      if (cflags.isSet(MNP_EVENT_RATE) && !okProbability(cflags.getValue(MNP_EVENT_RATE))) {
        cflags.setParseMessage("--" + MNP_EVENT_RATE + " must be a value in [0,1]");
        return false;
      }
      if (cflags.isSet(INS_EVENT_RATE) && !okProbability(cflags.getValue(INS_EVENT_RATE))) {
        cflags.setParseMessage("--" + INS_EVENT_RATE + " must be a value in [0,1]");
        return false;
      }
      if (cflags.isSet(DEL_EVENT_RATE) && !okProbability(cflags.getValue(DEL_EVENT_RATE))) {
        cflags.setParseMessage("--" + DEL_EVENT_RATE + " must be a value in [0,1]");
        return false;
      }

      if (!(cflags.isSet(COVERAGE) ^ cflags.isSet(READS))) {
        cflags.setParseMessage("Exactly one of --" + COVERAGE + " or --" + READS + " must be specified");
        return false;
      }
      if (cflags.isSet(READS)) {
        if ((Long) cflags.getValue(READS) <= 0) {
          cflags.setParseMessage("Number of reads should be greater than 0");
          return false;
        }
      } else if (cflags.isSet(COVERAGE)) {
        final double coverage = (Double) cflags.getValue(COVERAGE);
        if (coverage <= 0.0) {
          cflags.setParseMessage("Coverage should be positive");
          return false;
        } else if (coverage > 1000000.0) {
          cflags.setParseMessage("Coverage cannot be greater than 1000000.0");
          return false;
        }
      }
      if (!CommonFlags.validateOutputDirectory((File) cflags.getValue(OUTPUT))) {
        return false;
      }
      if (cflags.isSet(DISTRIBUTION) && cflags.isSet(TAXONOMY_DISTRIBUTION)) {
        cflags.setParseMessage("At most one of --" + DISTRIBUTION + " or --" + TAXONOMY_DISTRIBUTION + " should be set");
        return false;
      }
      if (cflags.isSet(DISTRIBUTION)) {
        final File distFile = (File) cflags.getValue(DISTRIBUTION);
        if (!distFile.exists() || distFile.isDirectory()) {
          cflags.setParseMessage("File: " + distFile + " does not exist");
          return false;
        }
      }
      if (!cflags.isSet(TAXONOMY_DISTRIBUTION) && (cflags.isSet(ABUNDANCE) || cflags.isSet(DNA_FRACTION))) {
        cflags.setParseMessage("--" + ABUNDANCE + " and --" + DNA_FRACTION + " are only applicable if using --" + TAXONOMY_DISTRIBUTION);
        return false;
      }
      if (cflags.isSet(TAXONOMY_DISTRIBUTION) && !(cflags.isSet(ABUNDANCE) || cflags.isSet(DNA_FRACTION))) {
        cflags.setParseMessage("either --" + ABUNDANCE + " or --" + DNA_FRACTION + " must be set if using --" + TAXONOMY_DISTRIBUTION);
        return false;

      }
      if (cflags.isSet(ABUNDANCE) && cflags.isSet(DNA_FRACTION)) {
        cflags.setParseMessage("At most one of --" + ABUNDANCE + " or --" + DNA_FRACTION + " should be set");
        return false;
      }
      if (cflags.isSet(TAXONOMY_DISTRIBUTION)) {
        final File distFile = (File) cflags.getValue(TAXONOMY_DISTRIBUTION);
        if (!distFile.exists() || distFile.isDirectory()) {
          cflags.setParseMessage("File: " + distFile + " does not exist");
          return false;
        }
      }

      if (cflags.isSet(QUAL_RANGE)) {
        final String range = (String) cflags.getValue(QUAL_RANGE);
        final String[] vals = range.split("-");
        if (vals.length != 2) {
          cflags.setParseMessage("Quality range is not of form qualmin-qualmax");
          return false;
        } else {
          try {
            final int l = Integer.parseInt(vals[0]);
            if ((l < 0) || (l > 63)) {
              cflags.setParseMessage("Minimum quality value must be between 0 and 63");
              return false;
            }
            final int u = Integer.parseInt(vals[1]);
            if ((u < 0) || (u > 63)) {
              cflags.setParseMessage("Maximum quality value must be between 0 and 63");
              return false;
            }
            if (l > u) {
              cflags.setParseMessage("Minimum quality value cannot be greater than maximum quality value");
              return false;
            }
          } catch (final NumberFormatException e) {
            cflags.setParseMessage("Quality range is not of form qualmin-qualmax");
            return false;
          }
        }
      }
      final File f = (File) cflags.getValue(INPUT);
      if (!f.exists()) {
        cflags.setParseMessage("The specified SDF, \"" + f.getPath() + "\", does not exist.");
        return false;
      }
      if (!f.isDirectory()) {
        cflags.setParseMessage("The specified file, \"" + f.getPath() + "\", is not an SDF.");
        return false;
      }
      if (cflags.isSet(TWIN_INPUT)) {
        final File tf = (File) cflags.getValue(TWIN_INPUT);
        if (!tf.exists()) {
          cflags.setParseMessage("The specified SDF, \"" + tf.getPath() + "\", does not exist.");
          return false;
        }
        if (!tf.isDirectory()) {
          cflags.setParseMessage("The specified file, \"" + tf.getPath() + "\", is not an SDF.");
          return false;
        }
        if (tf.equals(f)) {
          cflags.setParseMessage("The --" + TWIN_INPUT + " SDF cannot be the same as that given with --" + INPUT);
          return false;
        }
      }
      if ((Integer) cflags.getValue(MIN_FRAGMENT) > (Integer) cflags.getValue(MAX_FRAGMENT)) {
        cflags.setParseMessage("--" + MAX_FRAGMENT + " should not be smaller than --" + MIN_FRAGMENT);
        return false;
      }

      final File bed = (File) cflags.getValue(BED_FILE);
      if (cflags.isSet(BED_FILE)) {
        if (!bed.exists()) {
          cflags.setParseMessage("The --" + BED_FILE + " specified file doesn't exist: " + bed.getPath());
          return false;
        }
        if (cflags.isSet(COVERAGE)) {
          cflags.setParseMessage("--" + BED_FILE + " is incompatible with --" + COVERAGE);
          return false;
        }
      }

      return checkMachines(cflags);
    }

    protected boolean checkMachines(CFlags cflags) {
      final MachineType mt = MachineType.valueOf(cflags.getValue(MACHINE_TYPE).toString().toLowerCase(Locale.getDefault()));
      if (mt == MachineType.ILLUMINA_SE) {
        if (!cflags.checkRequired(READLENGTH)) {
          return false;
        }
        if ((Integer) cflags.getValue(READLENGTH) <= 1) {
          cflags.setParseMessage("Read length is too small");
          return false;
        }
        if ((Integer) cflags.getValue(READLENGTH) > (Integer) cflags.getValue(MIN_FRAGMENT)) {
          cflags.error("Read length is too large for selected fragment size");
          return false;
        }
        if (!cflags.checkBanned(LEFT_READLENGTH, RIGHT_READLENGTH, MIN_TOTAL_454_LENGTH, MAX_TOTAL_454_LENGTH, MIN_TOTAL_IONTORRENT_LENGTH, MAX_TOTAL_IONTORRENT_LENGTH)) {
          return false;
        }
      } else if (mt == MachineType.ILLUMINA_PE) {
        if (!cflags.checkRequired(LEFT_READLENGTH, RIGHT_READLENGTH)) {
          return false;
        }
        if (((Integer) cflags.getValue(LEFT_READLENGTH) <= 1)
            || ((Integer) cflags.getValue(RIGHT_READLENGTH) <= 1)) {
          cflags.error("Read length is too small");
          return false;
        }
        if (((Integer) cflags.getValue(LEFT_READLENGTH) > (Integer) cflags.getValue(MIN_FRAGMENT))
            || ((Integer) cflags.getValue(RIGHT_READLENGTH) > (Integer) cflags.getValue(MIN_FRAGMENT))) {
          cflags.error("Read length is too large for selected fragment size");
          return false;
        }
        if (!cflags.checkBanned(READLENGTH, MIN_TOTAL_454_LENGTH, MAX_TOTAL_454_LENGTH)) {
          return false;
        }

      } else if (mt == MachineType.COMPLETE_GENOMICS) {
        if (!cflags.checkBanned(READLENGTH, LEFT_READLENGTH, RIGHT_READLENGTH, MIN_TOTAL_454_LENGTH, MAX_TOTAL_454_LENGTH, MIN_TOTAL_IONTORRENT_LENGTH, MAX_TOTAL_IONTORRENT_LENGTH)) {
          return false;
        }

      } else if (mt == MachineType.FOURFIVEFOUR_PE || mt == MachineType.FOURFIVEFOUR_SE) {
        if (!cflags.checkRequired(MIN_TOTAL_454_LENGTH, MAX_TOTAL_454_LENGTH)) {
          return false;
        }
        if (!cflags.checkBanned(READLENGTH, LEFT_READLENGTH, RIGHT_READLENGTH, MIN_TOTAL_IONTORRENT_LENGTH, MAX_TOTAL_IONTORRENT_LENGTH)) {
          return false;
        }
        if ((Integer) cflags.getValue(MAX_TOTAL_454_LENGTH) > (Integer) cflags.getValue(MIN_FRAGMENT)) {
          cflags.error("Read length is too large for selected fragment size");
          return false;
        }
      } else if (mt == MachineType.IONTORRENT) {
        if (!cflags.checkRequired(MIN_TOTAL_IONTORRENT_LENGTH, MAX_TOTAL_IONTORRENT_LENGTH)) {
          return false;
        }
        if (!cflags.checkBanned(READLENGTH, LEFT_READLENGTH, RIGHT_READLENGTH, MIN_TOTAL_454_LENGTH, MAX_TOTAL_454_LENGTH)) {
          return false;
        }
        if ((Integer) cflags.getValue(MAX_TOTAL_IONTORRENT_LENGTH) > (Integer) cflags.getValue(MIN_FRAGMENT)) {
          cflags.error("Read length is too large for selected fragment size");
          return false;
        }
      } else {
        throw new IllegalArgumentException("Unhandled machine type: " + mt);
      }
      return true;
    }
  }

  private PortableRandom mRandom = null;
  private AbstractMachineErrorParams mPriors = null;

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return "generate simulated reads from a sequence";
  }

  @Override
  protected File outputDirectory() {
    File f = (File) mFlags.getValue(OUTPUT);
    if ("-".equals(f.getName())) {
      try {
        f = FileUtils.createTempDir("readsim", null);
        cleanDirectory();
      } catch (final IOException e) {
        throw new NoTalkbackSlimException("Could not create temporary directory " + e.getMessage());
      }
    }
    return f;
  }

  protected MachineType getMachineType() {
    return MachineType.valueOf(mFlags.getValue(MACHINE_TYPE).toString().toLowerCase(Locale.getDefault()));
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription("Generates reads from a reference genome.");
    mFlags.setCategories(UTILITY, new String[]{INPUT_OUTPUT, CAT_FRAGMENTS, CAT_ILLUMINA_PE, CAT_ILLUMINA_SE, CAT_454_PE, CAT_ION_SE, UTILITY});
    mFlags.registerExtendedHelp();
    mFlags.registerRequired('o', OUTPUT, File.class, "SDF", "name for reads output SDF").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('t', INPUT, File.class, "SDF", "SDF containing input genome").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional('T', TWIN_INPUT, File.class, "SDF", "SDF with second genome for simulating diploid genomes (deprecated)").setCategory(INPUT_OUTPUT);
    final Flag covFlag = mFlags.registerOptional('c', COVERAGE, Double.class, "float", "coverage, must be positive").setCategory(CAT_FRAGMENTS);
    final Flag nFlag = mFlags.registerOptional('n', READS, Long.class, "int", "number of reads to be generated").setCategory(CAT_FRAGMENTS);

    // Fragmenter
    mFlags.registerOptional('N', ALLOW_UNKNOWNS, "allow reads to be drawn from template fragments containing unknown nucleotides").setCategory(CAT_FRAGMENTS);
    mFlags.registerOptional('D', DISTRIBUTION, File.class, "file", "file containing probability distribution for sequence selection").setCategory(CAT_FRAGMENTS);
    mFlags.registerOptional(TAXONOMY_DISTRIBUTION, File.class, "file", "file containing probability distribution for sequence selection expressed as taxonomy id").setCategory(CAT_FRAGMENTS);
    mFlags.registerOptional(ABUNDANCE, "taxonomy distribution represents desired abundance").setCategory(CAT_FRAGMENTS);
    mFlags.registerOptional(DNA_FRACTION, "taxonomy distribution represents desired DNA fraction").setCategory(CAT_FRAGMENTS);
    mFlags.registerOptional(N_RATE, Double.class, "float", "rate that the machine will generate new unknowns in the read", 0.0).setCategory(CAT_FRAGMENTS);

    mFlags.registerOptional('s', SEED, Long.class, "int", "seed for random number generator").setCategory(UTILITY);
    mFlags.registerOptional(COMMENT, String.class, "string", "comment to include in the generated SDF").setCategory(UTILITY);
    SamCommandHelper.initSamRg(mFlags, "ILLUMINA", UTILITY);

    mFlags.addRequiredSet(covFlag);
    mFlags.addRequiredSet(nFlag);

    mFlags.registerOptional('q', QUAL_RANGE, String.class, "string", "set the range of base quality values permitted e.g.: 3-40. Default is fixed qualities corresponding to overall machine base error rate").setCategory(CommonFlagCategories.UTILITY);

    //reduce wastage
    mFlags.registerOptional(NO_NAMES, "do not create read names in result sdf").setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional(NO_QUAL, "do not create read qualities in result sdf").setCategory(CommonFlagCategories.UTILITY);

    // Override rate options
    mFlags.registerOptional(MNP_EVENT_RATE, Double.class, "float", "override the overall MNP event rate in the priors").setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional(INS_EVENT_RATE, Double.class, "float", "override the overall insertion event rate in the priors").setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional(DEL_EVENT_RATE, Double.class, "float", "override the overall deletion event rate in the priors").setCategory(CommonFlagCategories.UTILITY);
    // Limit to bed regions
    mFlags.registerOptional(BED_FILE, File.class, "FILE", "simulate exome capture by only generating reads that lie over the specified regions").setCategory(UTILITY);

    mFlags.registerOptional(PCR_DUP_RATE, Double.class, "float", "set the PCR duplication error rate", 0.0).setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional(CHIMERA_RATE, Double.class, "float", "set the chimeric fragment error rate", 0.0).setCategory(CommonFlagCategories.UTILITY);

    mFlags.setValidator(new ReadSimValidator());
    initMachineFlags();
  }
  protected void initMachineFlags() {
    initIlluminaFlags();
    init454Flags();
    initIonFlags();
    final Flag machType = mFlags.registerRequired(MACHINE_TYPE, String.class, "string", "select the sequencing technology to model").setCategory(INPUT_OUTPUT);
    machType.setParameterRange(new String[] {MachineType.ILLUMINA_SE.name(), MachineType.ILLUMINA_PE.name(), MachineType.COMPLETE_GENOMICS.name(), MachineType.FOURFIVEFOUR_PE.name(), MachineType.FOURFIVEFOUR_SE.name(), MachineType.IONTORRENT.name()});
    mFlags.registerOptional('E', MACHINE_ERROR_PRIORS, String.class, "string", "selects the sequencer machine error settings. One of [default, illumina, ls454_se, ls454_pe, complete, iontorrent]").setCategory(UTILITY);
  }

  protected void initIlluminaFlags() {
    // Illumina SE
    mFlags.registerOptional('r', READLENGTH, Integer.class, "int", "target read length, must be positive").setCategory(CAT_ILLUMINA_SE);
    mFlags.registerOptional('M', MAX_FRAGMENT, Integer.class, "int", "maximum fragment size", 250).setCategory(CAT_FRAGMENTS);
    mFlags.registerOptional('m', MIN_FRAGMENT, Integer.class, "int", "minimum fragment size", 200).setCategory(CAT_FRAGMENTS);
    // Illumina PE
    mFlags.registerOptional('L', LEFT_READLENGTH, Integer.class, "int", "target read length on the left side").setCategory(CAT_ILLUMINA_PE);
    mFlags.registerOptional('R', RIGHT_READLENGTH, Integer.class, "int", "target read length on the right side").setCategory(CAT_ILLUMINA_PE);
  }

  protected void init454Flags() {
    mFlags.registerOptional(MAX_TOTAL_454_LENGTH, Integer.class, "int", "maximum 454 read length (in paired end case the sum of the left and the right read lengths)").setCategory(CAT_454_PE);
    mFlags.registerOptional(MIN_TOTAL_454_LENGTH, Integer.class, "int", "minimum 454 read length (in paired end case the sum of the left and the right read lengths)").setCategory(CAT_454_PE);
  }

  protected void initIonFlags() {
    mFlags.registerOptional(MAX_TOTAL_IONTORRENT_LENGTH, Integer.class, "int", "maximum IonTorrent read length").setCategory(CAT_ION_SE);
    mFlags.registerOptional(MIN_TOTAL_IONTORRENT_LENGTH, Integer.class, "int", "minimum IonTorrent read length").setCategory(CAT_ION_SE);
  }

  private Machine createMachine() {
    final MachineType mt = getMachineType();
    final long seed = mRandom.nextLong();
    final Machine result;
    if (mt == MachineType.ILLUMINA_SE) {
      final IlluminaSingleEndMachine m = new IlluminaSingleEndMachine(mPriors, seed);
      m.setReadLength((Integer) mFlags.getValue(READLENGTH));
      result = m;
    } else if (mt == MachineType.ILLUMINA_PE) {
      final IlluminaPairedEndMachine m = new IlluminaPairedEndMachine(mPriors, seed);
      m.setLeftReadLength((Integer) mFlags.getValue(LEFT_READLENGTH));
      m.setRightReadLength((Integer) mFlags.getValue(RIGHT_READLENGTH));
      result = m;
    } else if (mt == MachineType.FOURFIVEFOUR_PE) {
      final FourFiveFourPairedEndMachine m = new FourFiveFourPairedEndMachine(mPriors, seed);
      m.setMinPairSize((Integer) mFlags.getValue(MIN_TOTAL_454_LENGTH));
      m.setMaxPairSize((Integer) mFlags.getValue(MAX_TOTAL_454_LENGTH));
      result = m;
    } else if (mt == MachineType.FOURFIVEFOUR_SE) {
      final FourFiveFourSingleEndMachine m = new FourFiveFourSingleEndMachine(mPriors, seed);
      m.setMinSize((Integer) mFlags.getValue(MIN_TOTAL_454_LENGTH));
      m.setMaxSize((Integer) mFlags.getValue(MAX_TOTAL_454_LENGTH));
      result = m;
    } else if (mt == MachineType.COMPLETE_GENOMICS) {
      result = new CompleteGenomicsMachine(mPriors, seed);
    } else if (mt == MachineType.IONTORRENT) {
      final IonTorrentSingleEndMachine m = new IonTorrentSingleEndMachine(mPriors, seed);
      m.setMinSize((Integer) mFlags.getValue(MIN_TOTAL_IONTORRENT_LENGTH));
      m.setMaxSize((Integer) mFlags.getValue(MAX_TOTAL_IONTORRENT_LENGTH));
      result = m;
    } else {
      throw new IllegalArgumentException("Unrecognized machine type: " + mt);
    }
    if (mFlags.isSet(QUAL_RANGE)) {
      final String range = (String) mFlags.getValue(QUAL_RANGE);
      final String[] vals = range.split("-");
      try {
        final int l = Integer.parseInt(vals[0]);
        final int u = Integer.parseInt(vals[1]);
        result.setQualRange((byte) l, (byte) u);
      } catch (final NumberFormatException e) {
        throw new NoTalkbackSlimException("Malformed quality range " + range);
      }
    }
    return result;
  }

  private ReadWriter createReadWriter(Machine m) throws IOException {
    final File f = (File) mFlags.getValue(OUTPUT);
    if (f.getName().endsWith(".fq")) {
      return new FastqReadWriter(f);
    } else if ("-".equals(f.getName())) {
      return new FastqReadWriter(System.out);
    } else {
      final SdfReadWriter rw = new SdfReadWriter(f, m.isPaired(), m.machineType(), !mFlags.isSet(NO_NAMES), !mFlags.isSet(NO_QUAL));
      rw.setComment((String) mFlags.getValue(COMMENT));
      if (mFlags.isSet(SamCommandHelper.SAM_RG)) {
        rw.setReadGroup(SamCommandHelper.validateAndCreateSamRG((String) mFlags.getValue(SamCommandHelper.SAM_RG), SamCommandHelper.ReadGroupStrictness.REQUIRED));
      }
      return rw;
    }
  }

  protected String getPriorsNameFlagValue() {
    if (mFlags.isSet(MACHINE_ERROR_PRIORS)) {
      return (String) mFlags.getValue(MACHINE_ERROR_PRIORS);
    }
    return null;
  }

  private AbstractMachineErrorParams createPriors() throws IOException {
    String priorsName = getPriorsNameFlagValue();
    if (priorsName == null) {
      final MachineType mt = getMachineType();
      priorsName = mt.priors();
    }
    try {
      final AbstractMachineErrorParams machineErrors = MachineErrorParams.builder().errors(priorsName).create();
      // Override rates if appropriate
      if (mFlags.isSet(MNP_EVENT_RATE) || mFlags.isSet(INS_EVENT_RATE) || mFlags.isSet(DEL_EVENT_RATE)) {
        final MachineErrorParamsBuilder mb = new MachineErrorParamsBuilder(machineErrors);
        if (mFlags.isSet(MNP_EVENT_RATE)) {
          mb.errorMnpEventRate((Double) mFlags.getValue(MNP_EVENT_RATE));
        }
        if (mFlags.isSet(INS_EVENT_RATE)) {
          mb.errorInsEventRate((Double) mFlags.getValue(INS_EVENT_RATE));
        }
        if (mFlags.isSet(DEL_EVENT_RATE)) {
          mb.errorDelEventRate((Double) mFlags.getValue(DEL_EVENT_RATE));
        }
        return mb.create();
      } else {
        return machineErrors;
      }
    } catch (final InvalidParamsException e) {
      return null;
    }
  }

  private Map<String, Double> createSelectionDistribution(final InputStream is) throws IOException {
    HashMap<String, Double> map = new HashMap<>();
    double sum = 0;
    try (BufferedReader r = new BufferedReader(new InputStreamReader(is))) {
      String line;
      while ((line = r.readLine()) != null) {
        if (line.length() > 0 && line.charAt(0) != '#') {
          final String[] parts = line.trim().split("\\s+");
          if (parts.length != 2) {
            throw new IOException("Malformed line: " + line);
          }
          try {
            final double p = Double.parseDouble(parts[0]);
            if (p < 0 || p > 1) {
              throw new IOException("Malformed line: " + line);
            }
            sum += p;
            final String species = parts[1].trim();
            if (map.containsKey(species)) {
              throw new IOException("Duplicated key: " + line);
            }
            map.put(species, p);
          } catch (final NumberFormatException e) {
            throw new IOException("Malformed line: " + line);
          }
        }
      }
    }
    if (Math.abs(sum - 1) > 0.00001) {
      Diagnostic.warning("Input distribution sums to: " + String.format("%1.5g", sum));
      final Map<String, Double> oldmap = map;
      map = new HashMap<>();
      for (Map.Entry<String, Double> entry : oldmap.entrySet()) {
        map.put(entry.getKey(), entry.getValue() / sum);
      }
    }
    return map;
  }

  void checkReadersForFragments(SequencesReader input, SequencesReader twinInput) {
    final int minFragment = (Integer) mFlags.getValue(MIN_FRAGMENT);
    final long minSequenceLength = twinInput == null ? input.minLength() : Math.min(input.minLength(), twinInput.minLength());
    final long maxSequenceLength = twinInput == null ? input.maxLength() : Math.max(input.maxLength(), twinInput.maxLength());
    if (minFragment > maxSequenceLength) {
      throw new NoTalkbackSlimException("All template sequences are too short for specified fragment lengths.");
    } else if (minFragment > minSequenceLength) {
      Diagnostic.warning(WarningType.INFO_WARNING, "The template contains some sequences that have length less than the minimum fragment length and these will not be present in the output.");
    }
  }
  @Override
  protected int mainExec(OutputStream out, LogStream log) throws IOException {
    final File input = (File) mFlags.getValue(INPUT);
    final File twinInput = (File) mFlags.getValue(TWIN_INPUT);
    try {
      final SequencesReader reader = SequencesReaderFactory.createMemorySequencesReaderCheckEmpty(input, true, false, LongRange.NONE);
      if (reader.numberSequences() > Integer.MAX_VALUE) {
        throw new NoTalkbackSlimException("Too many sequences");
      }
      final int numSeq = (int) reader.numberSequences();
      final double[] selectionProb;
      if (mFlags.isSet(DISTRIBUTION)) {
        Diagnostic.userLog("Using standard distribution");
        final double[] selectionDist = new double[numSeq];
        try (FileInputStream is = new FileInputStream((File) mFlags.getValue(DISTRIBUTION))) {
          final Map<String, Double> selectionMap = createSelectionDistribution(is);
          final PrereadNamesInterface names = reader.names();
          final int[] lengths = reader.sequenceLengths(0, numSeq);
          double sum = 0;
          for (int k = 0; k < numSeq; k++) {
            final Double p = selectionMap.get(names.name(k));
            if (p != null) {
              sum += p;
              selectionDist[k] = p * lengths[k];
            }
          }
          if (Math.abs(sum - 1) > 0.00001) {
            throw new NoTalkbackSlimException("Some sequences not seen in supplied template, sum:" + String.format("%1.5g", sum));
          }
          selectionProb = MathUtils.renormalize(selectionDist);
        }
        Diagnostic.userLog("Distribution complete");
      } else if (mFlags.isSet(TAXONOMY_DISTRIBUTION)) {
        Diagnostic.userLog("Using taxonomy distribution");
        final TaxonomyDistribution dist;
        try (final FileInputStream is = new FileInputStream((File) mFlags.getValue(TAXONOMY_DISTRIBUTION))) {
          dist = new TaxonomyDistribution(is, TaxonomyUtils.loadTaxonomyMapping(reader), reader, mFlags.isSet(DNA_FRACTION) ? TaxonomyDistribution.DistributionType.DNA_FRACTION : TaxonomyDistribution.DistributionType.ABUNDANCE);
        }
         selectionProb = dist.getDistribution();
        Diagnostic.userLog("Distribution complete");
      } else {
        selectionProb = null;
      }
      long totalResidues = reader.totalLength(); // subtract Ns
      try {
        final SequencesReader twinReader;
        if (mFlags.isSet(TWIN_INPUT)) {
          try {
            twinReader = SequencesReaderFactory.createMemorySequencesReaderCheckEmpty(twinInput, true, false, LongRange.NONE);
            if (twinReader.getSdfId().check(reader.getSdfId())) { //  twinReader.getSdfId() != 0 && twinReader.getSdfId() == reader.getSdfId()) {
              throw new NoTalkbackSlimException("The --" + TWIN_INPUT + " SDF cannot be the same as that given with --" + INPUT);
            }
            totalResidues += twinReader.totalLength();
            totalResidues /= 2; // When using a diploid genome we don't want to count both chromosome copies separately.
          } catch (final FileNotFoundException e) {
            throw new NoTalkbackSlimException(ErrorType.SDF_INDEX_NOT_VALID, twinInput.toString());
          }
        } else {
          twinReader = null;
        }
        checkReadersForFragments(reader, twinReader);
        try {
          if (reader.type() != SequenceType.DNA || twinReader != null && twinReader.type() != SequenceType.DNA) {
            throw new NoTalkbackSlimException("Input SDFs must be DNA");
          }
          if (mFlags.isSet(SEED)) {
            mRandom = new PortableRandom((Long) mFlags.getValue(SEED));
          } else {
            mRandom = new PortableRandom();
          }
          final long seed = mRandom.getSeed();
          mPriors = createPriors();
          if (mPriors == null) {
            mFlags.error(mFlags.getInvalidFlagMsg());
            cleanDirectory();
            return 1;
          }
          // Construct appropriate GenomeFragmenter / Machine / ReadWriter and validate
          final SequencesReader[] readers;
          final SequenceDistribution[] distributions;
          if (twinReader == null) {
            readers = new SequencesReader[] {reader};
            distributions = new SequenceDistribution[] {SequenceDistribution.createDistribution(reader, selectionProb)};
          } else {
            readers = new SequencesReader[] {reader, twinReader};
            distributions = new SequenceDistribution[] {SequenceDistribution.createDistribution(reader, selectionProb), SequenceDistribution.createDistribution(twinReader, selectionProb)};
          }

          final GenomeFragmenter gf = getGenomeFragmenter(readers, distributions);

          final Machine m;

          final Double pcrDupRate = (Double) mFlags.getValue(PCR_DUP_RATE);
          final Double chimeraRate = (Double) mFlags.getValue(CHIMERA_RATE);
          if (pcrDupRate > 0.0 || chimeraRate > 0.0) {
            m = new ErrorMachine(seed, createMachine(), pcrDupRate, chimeraRate);
          } else {
            m = createMachine();
          }

          Diagnostic.userLog("ReadSimParams" + LS + " input=" + input + LS + (twinInput == null ? "" : " diploid=" + twinInput + LS) + " machine=" + getMachineType() + LS + " output=" + outputDirectory() + LS + (mFlags.isSet(READS) ? " num-reads=" + mFlags.getValue(READS) + LS : "") + (mFlags.isSet(COVERAGE) ? " coverage=" + mFlags.getValue(COVERAGE) + LS : "") + (selectionProb == null ? "" : " distribution=" + Arrays.toString(selectionProb) + LS) + " allow-unknowns=" + mFlags.isSet(ALLOW_UNKNOWNS) + LS + " max-fragment=" + mFlags.getValue(MAX_FRAGMENT) + LS + " min-fragment=" + mFlags.getValue(MIN_FRAGMENT) + LS + " seed=" + seed + LS + LS + mPriors.toString() + LS);
          try (ReadWriter rw = getNFilter(createReadWriter(m))) {
            m.setReadWriter(rw);
            gf.setMachine(m);
            // Run generation
            if (mFlags.isSet(READS)) {
              fragmentByCount(gf, rw);
            } else {
              fragmentByCoverage(totalResidues, gf, m);
            }
            final double effectiveCoverage = (double) m.residues() / totalResidues;
            Diagnostic.info("Generated " + rw.readsWritten() + " reads, effective coverage " + Utils.realFormat(effectiveCoverage, 2));
            if (selectionProb != null) {
              FileUtils.stringToFile(gf.fractionStatistics(), new File(outputDirectory(), "fractions.tsv"));
            }
            //writeTemplateMappingFile(getTemplateMapping(reader, twinReader));
            Diagnostic.info(m.formatActionsHistogram());
          }
        } finally {
          if (twinReader != null) {
            twinReader.close();
          }
        }
      } catch (final FileNotFoundException e) {
        throw new NoTalkbackSlimException("Input SDF files are invalid.");
      } finally {
        reader.close();
      }
    } catch (final SecurityException e) {
      throw new IOException("Unable to create directory");
    } catch (final FileNotFoundException e) {
      throw new NoTalkbackSlimException(ErrorType.SDF_INDEX_NOT_VALID, input.toString());
    }
    return 0;
  }

  private GenomeFragmenter getGenomeFragmenter(SequencesReader[] readers, SequenceDistribution[] distributions) throws IOException {
    final GenomeFragmenter gf;
    if (mFlags.isSet(BED_FILE)) {
      final ReferenceRegions referenceRegions = BedUtils.regions((File) mFlags.getValue(BED_FILE));
      gf = new FilteringFragmenter(referenceRegions, mRandom.nextLong(), distributions, readers);
    } else {
      gf = new GenomeFragmenter(mRandom.nextLong(), distributions, readers);
    }
    gf.setMaxFragmentSize((Integer) mFlags.getValue(MAX_FRAGMENT));
    gf.setMinFragmentSize((Integer) mFlags.getValue(MIN_FRAGMENT));
    gf.allowNs(mFlags.isSet(ALLOW_UNKNOWNS));
    return gf;
  }

  private void fragmentByCoverage(long totalResidues, GenomeFragmenter gf, Machine m) throws IOException {
    final double targetCoverage = (Double) mFlags.getValue(COVERAGE);
    double coverage;
    long percentageDone = 0;
    do {
      gf.makeFragment();
      coverage = (double) m.residues() / totalResidues;
      final long percentage = (long) (coverage / targetCoverage * 100);
      if (percentage > percentageDone) {
        percentageDone = percentage;
        Diagnostic.progress(Math.min(100, percentageDone) + "% of reads generated");
      }
    } while (coverage < targetCoverage);
  }

  private void fragmentByCount(GenomeFragmenter gf, ReadWriter writer) throws IOException {
    final long targetReads = (Long) mFlags.getValue(READS);
    final long percentageIncrement = Math.max(targetReads / 100, 1);
    int lastIncrement = 0;
    int written;
    while ((written = writer.readsWritten()) < targetReads) {
      if (written - lastIncrement >= percentageIncrement) {
        Diagnostic.progress(((int) Math.min(100, written / (double) targetReads * 100.0)) + "% of reads generated");
        lastIncrement = written;
      }
      gf.makeFragment();
    }
  }

  private ReadWriter getNFilter(final ReadWriter internal) {
    final ReadWriter rw;
    final double nRate = (Double) mFlags.getValue(N_RATE);
    if (nRate > 0) {
      rw = new UnknownBaseReadWriter(internal, nRate, new PortableRandom(mRandom.getSeed()));
    } else {
      rw = internal;
    }
    return rw;
  }

  /**
   * Generate reads with specified coverage.
   * @param args arguments
   */
  public static void main(final String[] args) {
    new ReadSimCli().mainExit(args);
  }
}
