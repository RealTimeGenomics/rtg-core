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
package com.rtg.ngs;

import java.io.File;
import java.util.Locale;

import com.rtg.alignment.AlignerMode;
import com.rtg.alignment.EditDistanceFactory;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.HashingRegion;
import com.rtg.reader.FormatCli;
import com.rtg.reference.Sex;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.machine.MachineOrientation;

/**
 * Flags common between different map commands
 */
public final class MapFlags {

  /** Default word size */
  public static final int DEFAULT_WORD_SIZE = 22;
  /** score indel flag */
  public static final String XSCORE_INDEL = "Xscoreindel";
  /** flag for whether to compress hashes or not. */
  public static final String COMPRESS_HASHES_FLAG = "Xcompress-hashes";
  /** Max top results flag. */
  public static final String MAX_TOP_RESULTS_FLAG = "max-top-results";
  /** don't report short (un-indexable) reads */
  public static final String NO_SHORT = "Xdiscard-short-reads";
  /** don't report unmapped */
  public static final String NO_UNMAPPED = "no-unmapped";
  /** don't report unmated */
  public static final String NO_UNMATED = "no-unmated";
  /** Max allowed mismatches for unmated reads */
  public static final String UNMATED_MISMATCH_THRESHOLD = "max-unmated-mismatches";
  /** Max allowed mismatches for mated reads */
  public static final String MATED_MISMATCH_THRESHOLD = "max-mated-mismatches";
  /** Mask flag. */
  public static final String MASK_FLAG = "Xmask";
  /** Top N reads flag. */
  public static final String TOPN_RESULTS_FLAG = "Xmax-topn-results";
  /** Default top n value */
  public static final int DEFAULT_TOP_N = 10;
  /** <code>topequal</code> output filter */
  public static final String TOPEQUAL = "topequal";
  /** <code>topn</code> output filter */
  public static final String TOPN = "topn";
  /** max matches allowed during alignment flag */
  public static final String MAX_ALIGNMENT_MISMATCHES = "max-mismatches";
  /** flag for setting the step size */
  public static final String STEP_FLAG = "step";
  /** Specify whether temporary files should be compressed */
  public static final String TEMP_FILES_COMPRESSED = "Xtemp-files-gzipped";
  /** indel length flag */
  public static final String INDEL_LENGTH_FLAG = "indel-length";
  /** substitution flag */
  public static final String SUBSTITUTIONS_FLAG = "substitutions";
  /** word size flag */
  public static final String WORDSIZE_FLAG = "word";
  /** Blacklist threshold flag */
  public static final String BLACKLIST_THRESHOLD = "blacklist-threshold";
  /** Maximum repeat frequency */
  public static final String MAX_REPEAT_FREQUENCY_FLAG = "Xmax-repeat-freq";
  /** Maximum repeat frequency */
  public static final String MIN_REPEAT_FREQUENCY_FLAG = "Xmin-repeat-freq";
  /** Force mapping commands to use long read code path instead of short read */
  public static final String FORCE_LONG_FLAG = "Xforce-long";
  /** Parallel unmated processing */
  public static final String PARALLEL_UNMATED_PROCESSING_FLAG = "Xparallel-unmated-processing";
  private static final int MAX_WORD_SIZE = 32; //Integer.valueOf(System.getProperty("rtg.max-word", "32"));
  private static final int MAX_INSERT_SIZE = 1000000000;
  private static final int BITS_FOR_SCORE1 = 12;
  /** Maximum value of score that can be stored in top n data structure. */
  public static final int MAX_SCORE = (1 << BITS_FOR_SCORE1) - 1;

  private MapFlags() { }

  static final String[] FORMAT_OPTIONS = {FormatCli.SDF_FORMAT, FormatCli.FASTA_FORMAT, FormatCli.FASTQ_FORMAT, FormatCli.SAM_SE_FORMAT, FormatCli.SAM_PE_FORMAT};


  //private static final String INSERT_SIZE_FLAG = "insert-size";
  static final String THREAD_MULTIPLIER = "Xthread-multiplier";
  static final String X_LONG_READ = "Xlong-read";
  static final String OUTPUT_UNFILTERED = "all-hits";
  static final String OUTPUT_NULLFILTERED = "Xnull-filter";
  static final String READ_FREQUENCY_FLAG = "Xread-freq";
  static final String MIN_HITS_FLAG = "Xmin-hits";
  static final String NO_INMEMORY_TEMPLATE = "Xno-inmemory-template";
  static final String OUTPUT_READ_NAMES_FLAG = "read-names";
  static final String LEGACY_CIGARS = "legacy-cigars";
  static final String BAM_FLAG = "bam";
  /** unify alignment outputs flag */
  public static final String DONT_UNIFY_FLAG = "no-merge";
  /** sam output flag */
  public static final String SAM_FLAG = "sam";
  /** Sex specification flag. */
  public static final String SEX_FLAG = "sex";
  /** Disable calibration */
  public static final String NO_CALIBRATION = "no-calibration";
  static final String NO_SVPREP = "no-svprep";

  static final String N_AS_MISMATCH = "penalize-unknowns";

  static final String GAP_OPEN_PENALTY_FLAG = "gap-open-penalty";
  static final String GAP_EXTEND_PENALTY_FLAG = "gap-extend-penalty";
  static final String MISMATCH_PENALTY_FLAG = "mismatch-penalty";
  static final String UNKNOWNS_PENALTY_FLAG = "unknowns-penalty";
  static final String SOFT_CLIP_DISTANCE_FLAG = "soft-clip-distance";
  /** Controls length/number of INDELS that can be handled during alignment. */
  static final String ALIGNER_BAND_WIDTH_FACTOR_FLAG = "aligner-band-width";

  static final String ALIGNER_MODE_FLAG = "aligner-mode";
  static final String SINGLE_INDEL_PENALTIES_FLAG = "Xsingle-indel-penalties";

  static void initMapFlags(CFlags flags) {
    flags.registerOptional('e', MAX_ALIGNMENT_MISMATCHES, IntegerOrPercentage.class, CommonFlags.INT, "maximum mismatches for mappings in single-end mode (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf(NgsFilterParams.MAX_MATED_MISMATCH_THRESHOLD)).setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional('E', UNMATED_MISMATCH_THRESHOLD, IntegerOrPercentage.class, CommonFlags.INT, "maximum mismatches for mappings of unmated results (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf(NgsFilterParams.MAX_UNMATED_MISMATCH_THRESHOLD)).setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(MATED_MISMATCH_THRESHOLD, IntegerOrPercentage.class, CommonFlags.INT, "maximum mismatches for mappings across mated results, alias for --" + MAX_ALIGNMENT_MISMATCHES + " (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf(NgsFilterParams.MAX_MATED_MISMATCH_THRESHOLD)).setCategory(CommonFlagCategories.REPORTING);

    flags.registerOptional(CommonFlags.TEMP_DIR, File.class, CommonFlags.DIR, "directory used for temporary files (Defaults to output directory)").setCategory(CommonFlagCategories.UTILITY);

    flags.registerOptional(LEGACY_CIGARS, "use legacy cigars in output").setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(OUTPUT_READ_NAMES_FLAG, "use read name in output instead of read id (Uses more RAM)").setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(ALIGNER_MODE_FLAG, AlignerMode.class, CommonFlags.STRING, "pick the aligner to use", AlignerMode.AUTO).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(SINGLE_INDEL_PENALTIES_FLAG, String.class, "STRING|FILE", "single indel penalty file", EditDistanceFactory.DEFAULT_SINGLE_INDEL_TABLE).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);

    //--X flags


    flags.registerOptional(THREAD_MULTIPLIER, Integer.class, CommonFlags.INT, "number of work chunks per thread", HashingRegion.DEFAULT_THREAD_MULTIPLIER).setCategory(CommonFlagCategories.UTILITY);
    //flags.registerOptional(INSERT_SIZE_FLAG, Integer.class, INT, "expected insert size for pairs");
    flags.registerOptional(TEMP_FILES_COMPRESSED, Boolean.class, "BOOL", "gzip temporary SAM files", Boolean.TRUE).setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(NO_INMEMORY_TEMPLATE, "do not load the template in memory").setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(FORCE_LONG_FLAG, "force the use of long read mode").setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(SEX_FLAG, Sex.class, "sex", "sex of sample", null).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(CommonFlags.PEDIGREE_FLAG, File.class, CommonFlags.FILE, "genome relationships pedigree containing sex of sample").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);

    CommonFlags.initReadRange(flags);
  }

  /**
   * Adds map related flags for aligner penalties
   * @param flags flags to set
   */
  public static void initAlignerPenaltyFlags(CFlags flags) {
    flags.registerOptional(MapFlags.GAP_OPEN_PENALTY_FLAG, Integer.class, CommonFlags.INT, "penalty for a gap open during alignment", EditDistanceFactory.DEFAULT_GAP_OPEN_PENALTY).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(MapFlags.GAP_EXTEND_PENALTY_FLAG, Integer.class, CommonFlags.INT, "penalty for a gap extension during alignment", EditDistanceFactory.DEFAULT_GAP_EXTEND_PENALTY).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(MapFlags.MISMATCH_PENALTY_FLAG, Integer.class, CommonFlags.INT, "penalty for a mismatch during alignment", EditDistanceFactory.DEFAULT_SUBSTITUTION_PENALTY).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(MapFlags.UNKNOWNS_PENALTY_FLAG, Integer.class, CommonFlags.INT, "penalty for unknown nucleotides during alignment", EditDistanceFactory.DEFAULT_UNKNOWNS_PENALTY).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(MapFlags.SOFT_CLIP_DISTANCE_FLAG, Integer.class, CommonFlags.INT, "soft clip alignments if indels occur INT bp from either end", 5).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);

    flags.register(new Flag<>(null, MapFlags.ALIGNER_BAND_WIDTH_FACTOR_FLAG, "aligner indel band width scaling factor, fraction of read length allowed as an indel", 0, 1, Double.class, CommonFlags.FLOAT, 0.5, CommonFlagCategories.SENSITIVITY_TUNING));
  }

  /**
   * adds map related input and output flags
   * @param flags flags to set
   */
  public static void initInputOutputFlags(CFlags flags) {
    final Flag<File> input = flags.registerOptional('i', CommonFlags.READS_FLAG, File.class, CommonFlags.SDF_OR_FILE, "input read set").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    CommonFlags.initOutputDirFlag(flags);
    initTemplateFlag(flags);

    initInputFormatFlags(flags);
    flags.addRequiredSet(input);
  }

  /**
   * Initialise mapping IO flags (input, output, template)
   * @param flags shared flags
   */
  public static void initMapIOFlags(CFlags flags) {
    flags.registerRequired('i', CommonFlags.READS_FLAG, File.class, CommonFlags.SDF, "SDF containing reads to map").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    initTemplateFlag(flags);
    CommonFlags.initOutputDirFlag(flags);
  }

  /**
   * Initialise the input format flags for mapping commands
   * @param flags shared flags
   */
  public static void initInputFormatFlags(CFlags flags) {
    //mapx already has a -f flag
    final Flag<String> formatFlag = flags.registerOptional('F', FormatCli.FORMAT_FLAG, String.class, "FORMAT", "input format for reads", FormatCli.SDF_FORMAT).setCategory(CommonFlagCategories.INPUT_OUTPUT);
    formatFlag.setParameterRange(FORMAT_OPTIONS);
    CommonFlags.initQualityFormatFlag(flags);
  }

  static void initPairedEndFormatFlags(CFlags flags)  {
    //IO format flags specific to paired end
    final Flag<File> inputL = flags.registerOptional('l', FormatCli.LEFT_FILE_FLAG, File.class, CommonFlags.FILE, "left input file for FASTA/FASTQ paired end reads").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    final Flag<File> inputR = flags.registerOptional('r', FormatCli.RIGHT_FILE_FLAG, File.class, CommonFlags.FILE, "right input file for FASTA/FASTQ paired end reads").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.addRequiredSet(inputL, inputR);
  }

  static void initReadFreqFlag(CFlags flags, int def) {
    flags.registerOptional(READ_FREQUENCY_FLAG, Integer.class, CommonFlags.INT, "maximum number of hits before a read is discarded", def).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
  }

  /**
   * Initialise the SAM output flag
   * @param flags shared flags
   */
  public static void initSamOutputFlag(CFlags flags) {
    flags.registerOptional(SAM_FLAG, "output the alignment files in SAM format").setCategory(CommonFlagCategories.INPUT_OUTPUT);
  }

  static void initNoCalibrationFlag(CFlags flags) {
    flags.registerOptional(NO_CALIBRATION, "do not produce calibration files").setCategory(CommonFlagCategories.UTILITY);
  }

  static void initSvPrepFlag(CFlags flags) {
    flags.registerOptional(NO_SVPREP, "do not perform structural variant processing").setCategory(CommonFlagCategories.UTILITY);
  }

  /**
   * Initialise the BAM output flag
   * @param flags shared flags
   * @param desc the flag description
   */
  static void initBamOutputFlag(CFlags flags, String desc) {
    flags.registerOptional(BAM_FLAG, desc).setCategory(CommonFlagCategories.INPUT_OUTPUT);
  }

  /**
   * Initialise the flag which disables unified alignment output
   * @param flags shared flags
   */
  static void initDontUnifyFlag(CFlags flags) {
    flags.registerOptional(DONT_UNIFY_FLAG, "output mated/unmated/unmapped alignments into separate SAM/BAM files").setCategory(CommonFlagCategories.UTILITY);
  }

  /**
   * Validate paired end flags for formatting options
   * @param flags shared flags
   * @return true if valid
   */
  public static boolean validatePairedEndFormatFlags(CFlags flags) {
    if (flags.isSet(CommonFlags.READS_FLAG)) {
      if (flags.isSet(FormatCli.LEFT_FILE_FLAG) || flags.isSet(FormatCli.RIGHT_FILE_FLAG)) {
        flags.setParseMessage("Can only specify --" + FormatCli.LEFT_FILE_FLAG + " and --" + FormatCli.RIGHT_FILE_FLAG + " when not using --" + CommonFlags.READS_FLAG);
        return false;
      }
    } else if ((flags.isSet(FormatCli.LEFT_FILE_FLAG) && !flags.isSet(FormatCli.RIGHT_FILE_FLAG))
        || (flags.isSet(FormatCli.RIGHT_FILE_FLAG) && !flags.isSet(FormatCli.LEFT_FILE_FLAG))) {
      flags.setParseMessage("Must set both --" + FormatCli.LEFT_FILE_FLAG + " and --" + FormatCli.RIGHT_FILE_FLAG);
      return false;
    } else {
      final String format = flags.isSet(FormatCli.FORMAT_FLAG) ? flags.getValue(FormatCli.FORMAT_FLAG).toString().toLowerCase(Locale.getDefault()) : FormatCli.SDF_FORMAT;
      if (FormatCli.SAM_SE_FORMAT.equals(format) || FormatCli.SAM_PE_FORMAT.equals(format)) {
        flags.setParseMessage("Do not use left and right flags when using SAM single or paired end input.");
        return false;
      }
    }
    return true;
  }

  static boolean validateMapParams(CFlags flags) {
    final int a = (Integer) flags.getValue(SUBSTITUTIONS_FLAG);
    final int b = (Integer) flags.getValue(CommonFlags.INDELS_FLAG);
    final int c = (Integer) flags.getValue(INDEL_LENGTH_FLAG);
    if (flags.isSet(WORDSIZE_FLAG)) {
      if (!flags.checkInRange(WORDSIZE_FLAG, 1, MAX_WORD_SIZE)) {
        return false;
      }
    }
    if (!validateStepAndWordSize(flags)) {
      return false;
    }
    if (a < 0) {
      Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "-a", a + "", "0");
      return false;
    }
    if (b < 0) {
      Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "-b", b + "", "0");
      return false;
    }
    if (c <= 0) {
      Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "-c", c + "", "1");
      return false;
    }
    if (!checkPercentRepeatFrequency(flags)) {
      return false;
    }
    if (!checkBlacklistThreshold(flags)) {
      return false;
    }

    if (!CommonFlags.validateThreads(flags)) {
      return false;
    }

    if (!validateMinMaxFragmentSize(flags)) {
      return false;
    }
    if (flags.isSet(MATED_MISMATCH_THRESHOLD) && flags.isSet(MAX_ALIGNMENT_MISMATCHES)) {
      flags.setParseMessage("--" + MATED_MISMATCH_THRESHOLD + " is an alias for --" + MAX_ALIGNMENT_MISMATCHES + ", please only set one");
      return false;
    }
    final IntegerOrPercentage maxMatedScore = (IntegerOrPercentage) flags.getValue(MATED_MISMATCH_THRESHOLD);
    if (maxMatedScore.getValue(100) < 0) {
      Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + MATED_MISMATCH_THRESHOLD, maxMatedScore + "", "0");
      return false;
    }
    final IntegerOrPercentage maxAlignmentScore = (IntegerOrPercentage) flags.getValue(MAX_ALIGNMENT_MISMATCHES);
    if (maxAlignmentScore.getValue(100) < 0) {
      Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + MAX_ALIGNMENT_MISMATCHES, maxAlignmentScore + "", "0");
      return false;
    }
    final IntegerOrPercentage maxUnMatedScore = (IntegerOrPercentage) flags.getValue(UNMATED_MISMATCH_THRESHOLD);
    if (maxUnMatedScore.getValue(100) < 0) {
      Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + UNMATED_MISMATCH_THRESHOLD, maxUnMatedScore + "", "0");
      return false;
    }

    if (!validatePenaltyFlags(flags)) {
      return false;
    }

    return true;
  }

  /**
   * Validate fragment size limiting flags
   * @param flags the flags to check
   * @return true if flags passed check
   */
  public static boolean validateMinMaxFragmentSize(CFlags flags) {
    if (!flags.checkInRange(CommonFlags.MAX_FRAGMENT_SIZE, 1, MAX_INSERT_SIZE)) {
      return false;
    }
    if (!flags.checkInRange(CommonFlags.MIN_FRAGMENT_SIZE, 0, MAX_INSERT_SIZE)) {
      return false;
    }
    if ((Integer) flags.getValue(CommonFlags.MIN_FRAGMENT_SIZE) > (int) (Integer) flags.getValue(CommonFlags.MAX_FRAGMENT_SIZE)) {
      flags.setParseMessage("--" + CommonFlags.MIN_FRAGMENT_SIZE + " must be less than --" + CommonFlags.MAX_FRAGMENT_SIZE);
      return false;
    }
    return true;
  }

  static boolean validateMapInputOutputParams(CFlags flags) {
    final String format = flags.isSet(FormatCli.FORMAT_FLAG) ? flags.getValue(FormatCli.FORMAT_FLAG).toString().toLowerCase(Locale.getDefault()) : FormatCli.SDF_FORMAT;
    final boolean sdf = FormatCli.SDF_FORMAT.equals(format);
    if (sdf) {
      if (!flags.isSet(CommonFlags.READS_FLAG)) {
        flags.setParseMessage("Must set --" + CommonFlags.READS_FLAG + " when --" + FormatCli.FORMAT_FLAG + " is set to " + FormatCli.SDF_FORMAT);
        return false;
      }
    } else {
      if (!flags.isSet(CommonFlags.READS_FLAG) && !flags.isSet(FormatCli.LEFT_FILE_FLAG) && !flags.isSet(FormatCli.RIGHT_FILE_FLAG)) {
        flags.setParseMessage("Must set reads flag. Either --" + CommonFlags.READS_FLAG + " or --" + FormatCli.LEFT_FILE_FLAG + " and --" + FormatCli.RIGHT_FILE_FLAG);
        return false;
      }
    }
    if (!FormatCli.validateQualityFormatFlags(flags, format)) {
      return false;
    }
    if (!validatePairedEndFormatFlags(flags)) {
      return false;
    }
    if (!CommonFlags.validateReads(flags, sdf)) {
      return false;
    }

    if (!CommonFlags.validateOutputDirectory(flags)) {
      return false;
    }
    if (!CommonFlags.validateTemplate(flags)) {
      return false;
    }
    //final File input = (File) flags.getValue(CommonFlags.READS_FLAG);
    final boolean isPaired = MapParamsHelper.isPaired(flags); //ReaderUtils.isPairedEndDirectory(input);
    if (flags.isSet(CommonFlags.MIN_FRAGMENT_SIZE) && !isPaired) {
      flags.setParseMessage("--" + CommonFlags.MIN_FRAGMENT_SIZE + " requires paired end input");
      return false;
    }
    if (flags.isSet(CommonFlags.MAX_FRAGMENT_SIZE) && !isPaired) {
      flags.setParseMessage("--" + CommonFlags.MAX_FRAGMENT_SIZE + " requires paired end input");
      return false;
    }
    if (flags.isSet(MATED_MISMATCH_THRESHOLD) && !isPaired) {
      flags.setParseMessage("--" + MATED_MISMATCH_THRESHOLD + " requires paired end input");
      return false;
    }
    if (flags.isSet(UNMATED_MISMATCH_THRESHOLD) && !isPaired) {
      flags.setParseMessage("--" + UNMATED_MISMATCH_THRESHOLD + " requires paired end input");
      return false;
    }

    if (flags.isSet(BAM_FLAG) && flags.isSet(SAM_FLAG)) {
      flags.setParseMessage("Can not specify both --" + SAM_FLAG + " and --" + BAM_FLAG);
      return false;
    }
    if (!flags.isSet(SAM_FLAG) && flags.isSet(CommonFlags.NO_GZIP)) {
      flags.setParseMessage("Can only specify --" + CommonFlags.NO_GZIP + " when using SAM output");
      return false;
    }
    return true;
  }

  static boolean validateSexTemplateReference(CFlags flags) {
    return CommonFlags.validateSexTemplateReference(flags, MapFlags.SEX_FLAG, CommonFlags.PEDIGREE_FLAG, CommonFlags.TEMPLATE_FLAG);
  }

  static boolean validatePenaltyFlags(CFlags flags) {
    if (!flags.checkInRange(MISMATCH_PENALTY_FLAG, 1, 100)) {
      return false;
    }
    if (!flags.checkInRange(UNKNOWNS_PENALTY_FLAG, 0, 100)) {
      return false;
    }
    if (!flags.checkInRange(GAP_OPEN_PENALTY_FLAG, 1, 100)) {
      return false;
    }
    if (!flags.checkInRange(GAP_EXTEND_PENALTY_FLAG, 1, 100)) {
      return false;
    }
    if (!flags.checkInRange(SOFT_CLIP_DISTANCE_FLAG, 0, 100)) {
      return false;
    }
    if (flags.isSet(ALIGNER_BAND_WIDTH_FACTOR_FLAG) && !flags.checkInRange(ALIGNER_BAND_WIDTH_FACTOR_FLAG, (double) 0.0f, (double) 1.0f)) {
      return false;
    }
    if (!(AlignerMode.GENERAL == flags.getValue(ALIGNER_MODE_FLAG)) && (flags.isSet(MISMATCH_PENALTY_FLAG) || flags.isSet(UNKNOWNS_PENALTY_FLAG) || flags.isSet(GAP_OPEN_PENALTY_FLAG) || flags.isSet(GAP_EXTEND_PENALTY_FLAG))) {
      flags.setParseMessage("Penalty flags are only valid if --" + ALIGNER_MODE_FLAG + " is 'general'");
      return false;

    }

    return true;
  }

  /**
   * Initialise shared flags
   * @param flags shared flags
   */
  public static void initSharedFlags(CFlags flags) {
    initMaskFlags(flags);
    initSharedFlagsOnly(flags);
    initStepSize(flags);
  }

  /**
   * Initialise shared flags without mask flags
   * @param flags shared flags
   */
  public static void initSharedFlagsOnly(CFlags flags) {
    initSharedFlagsOnly(flags, IntegerOrPercentage.valueOf("90%"), 1, 1000);
  }

  /**
   * Initialise shared flags without mask flags
   * @param flags shared flags
   * @param repeat value to set default repeat frequency
   * @param minRepeat value to set default minimum repeat frequency
   * @param maxRepeat value to set default max repeat frequency
   */
  public static void initSharedFlagsOnly(CFlags flags, IntegerOrPercentage repeat, int minRepeat, int maxRepeat) {
    flags.registerOptional(CommonFlags.REPEAT_FREQUENCY_FLAG, IntegerOrPercentage.class, CommonFlags.INT, "maximum repeat frequency", repeat).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(BLACKLIST_THRESHOLD, Integer.class, CommonFlags.INT, "filter k-mers that occur more than this many times in the reference using a blacklist").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(MAX_REPEAT_FREQUENCY_FLAG, Integer.class, CommonFlags.INT, "upper limit for repeat frequency when using the proportional repeat frequency setting", maxRepeat).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(MIN_REPEAT_FREQUENCY_FLAG, Integer.class, CommonFlags.INT, "lower limit for repeat frequency when using the proportional repeat frequency setting", minRepeat).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    CommonFlags.initNoGzip(flags);
    flags.registerOptional(PARALLEL_UNMATED_PROCESSING_FLAG, Boolean.class, CommonFlags.BOOL, "run unmated processing in parallel with mated processing", Boolean.FALSE).setCategory(CommonFlagCategories.UTILITY);
    CommonFlags.initThreadsFlag(flags);
  }

  /**
   * Init word size flag with given description.
   * @param flags flags to add to
   * @param desc description to register
   */
  public static void initWordSize(CFlags flags, String desc) {
    flags.registerOptional('w', WORDSIZE_FLAG, Integer.class, CommonFlags.INT, desc).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
  }

  /**
   * Init step size flag
   * @param flags flags to add to.
   */
  public static void initStepSize(CFlags flags) {
    initStepSize(flags, "step size (Default is word size)");
  }

  /**
   * Init step size flag
   * @param flags flags to add to.
   * @param desc description for flag
   */
  public static void initStepSize(CFlags flags, String desc) {
    flags.registerOptional('s', STEP_FLAG, Integer.class, CommonFlags.INT, desc).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
  }

  /**
   * Adds mask flags -w -a -b -c to the given <code>CFlags</code> object
   * @param flags flags to add to
   */
  public static void initMaskFlags(CFlags flags) {
    initWordSize(flags, "word size (Default is " + DEFAULT_WORD_SIZE + ", or read length / 2, whichever is smaller)");
    initMaskFlagsOnly(flags);
  }

  /**
   * Adds mask flags -a -b -c to the given <code>CFlags</code> object
   * @param flags flags to add to
   */
  public static void initMaskFlagsOnly(CFlags flags) {
    flags.registerOptional('a', SUBSTITUTIONS_FLAG, Integer.class, CommonFlags.INT, "guaranteed minimum number of substitutions which will be detected", 1).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional('b', CommonFlags.INDELS_FLAG, Integer.class, CommonFlags.INT, "guaranteed minimum number of indels which will be detected", 1).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional('c', INDEL_LENGTH_FLAG, Integer.class, CommonFlags.INT, "guaranteed number of positions that will be detected in a single indel", 1).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
  }

  /**
   * Initialise flags for paired end specific details (min and max insert sizes)
   * @param flags shared flags
   */
  public static void initFragmentSizeFlags(final CFlags flags) {
    flags.registerOptional('M', CommonFlags.MAX_FRAGMENT_SIZE, Integer.class, CommonFlags.INT, "maximum permitted fragment size when mating paired reads", 1000).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional('m', CommonFlags.MIN_FRAGMENT_SIZE, Integer.class, CommonFlags.INT, "minimum permitted fragment size when mating paired reads", 0).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
  }

  /**
   * Initialise flags for paired end specific details (min and max insert sizes)
   * @param flags shared flags
   */
  public static void initPairedEndFlags(final CFlags flags) {
    flags.registerOptional('d', CommonFlags.PAIR_ORIENTATION_FLAG, MachineOrientation.class, CommonFlags.STRING, "orientation for proper pairs", MachineOrientation.ANY).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    initFragmentSizeFlags(flags);
  }

  /**
   * This method checks that the repeat frequency flag is within appropriate boundaries.
   * @param flags the flags to check
   * @return <code>true</code> if all okay <code>false</code> otherwise
   */
  public static boolean checkRepeatFrequency(CFlags flags) {
    if (flags.isSet(CommonFlags.REPEAT_FREQUENCY_FLAG)) {
      if (!flags.checkInRange(CommonFlags.REPEAT_FREQUENCY_FLAG, 1, 100000)) {
        return false;
      }
    }
    return true;
  }

  /**
   * This method checks that the repeat frequency flag is within appropriate boundaries.
   * @param flags the flags to check
   * @return <code>true</code> if all okay <code>false</code> otherwise
   */
  public static boolean checkPercentRepeatFrequency(CFlags flags) {
    if (flags.isSet(CommonFlags.REPEAT_FREQUENCY_FLAG)) {
      final IntegerOrPercentage repeat = (IntegerOrPercentage) flags.getValue(CommonFlags.REPEAT_FREQUENCY_FLAG);
      if (repeat.isPercentage()) {
        if (repeat.getValue(100) < 0) {
          Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + CommonFlags.REPEAT_FREQUENCY_FLAG, repeat + "", "0%");
          return false;
        } else if (repeat.getValue(100) > 100) {
          Diagnostic.error(ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, "--" + CommonFlags.REPEAT_FREQUENCY_FLAG, repeat + "", "100%");
          return false;
        }
      } else {
        if (repeat.getValue(100) < 1) {
          Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + CommonFlags.REPEAT_FREQUENCY_FLAG, repeat + "", 1 + "");
          return false;
        }
        if (repeat.getValue(100) > 100000) {
          Diagnostic.error(ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, "--" + CommonFlags.REPEAT_FREQUENCY_FLAG, repeat + "", 100000 + "");
          return false;
        }
      }
    }
    if (!flags.checkInRange(MAX_REPEAT_FREQUENCY_FLAG, -1, Integer.MAX_VALUE)) {
      return false;
    }
    if (!flags.checkInRange(MIN_REPEAT_FREQUENCY_FLAG, -1, Integer.MAX_VALUE)) {
      return false;
    }
    return true;
  }

  /**
   * Checks related to blacklist threshold
   * @param flags the flags to check
   * @return <code>true</code> if all okay <code>false</code> otherwise
   */
  public static boolean checkBlacklistThreshold(CFlags flags) {
    if (flags.isSet(BLACKLIST_THRESHOLD)) {
      final int val = (Integer) flags.getValue(BLACKLIST_THRESHOLD);
      if (val < 0) {
        Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + CommonFlags.REPEAT_FREQUENCY_FLAG, val + "", "0");
        return false;
      }
    }
    return true;
  }

  /**
   * Initialise the template flag
   * @param flags shared flags
   */
  public static void initTemplateFlag(CFlags flags) {
    flags.registerRequired('t', CommonFlags.TEMPLATE_FLAG, File.class, CommonFlags.SDF, "SDF containing template to map against").setCategory(CommonFlagCategories.INPUT_OUTPUT);
  }

  /**
   * Return the word size based on word size flag or read length
   * @param flags the flags
   * @param readLength the read length
   * @param defWordSize the default max word size
   * @return the word size
   */
  public static int getWordSize(CFlags flags, int readLength, int defWordSize) {
    if (flags.isSet(WORDSIZE_FLAG)) {
      return (Integer) flags.getValue(WORDSIZE_FLAG);
    } else {
      return Math.min(defWordSize, readLength / 2);
    }
  }

  /**
   * Validate step and word size flags
   * @param flags shared flags
   * @return true if valid, otherwise false
   */
  public static boolean validateStepAndWordSize(CFlags flags) {
    final Integer wordSize = flags.isSet(WORDSIZE_FLAG) ? (Integer) flags.getValue(WORDSIZE_FLAG) : null;
    if (wordSize != null && wordSize < 1) {
      Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + WORDSIZE_FLAG, Integer.toString(wordSize), "1");
      return false;
    }
    if (flags.isSet(STEP_FLAG)) {
      final int step = (Integer) flags.getValue(STEP_FLAG);
      if (step < 1) {
        Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + STEP_FLAG, Integer.toString(step), "1");
        return false;
      }
    }
    return true;
  }

  /**
   * validate mask parameter
   * @param readLength read length
   * @param w word
   * @param subs substitution
   * @param i indel
   * @param l indel length
   * @throws InvalidParamsException if parameter is not valid
   */
  public static void validateMaskParams(final int readLength, final int w, final int subs, final int i, final int l) {
    if (w <= 0) {
      throw new InvalidParamsException(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "-w", w + "", "1");
    }
    if (subs < 0) {
      throw new InvalidParamsException(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "-a", subs + "", "0");
    }
    if (i < 0) {
      throw new InvalidParamsException(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "-b", i + "", "0");
    }
    if (l <= 0) {
      throw new InvalidParamsException(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "-c", l + "", "1");
    }
    if (readLength < w) {
      throw new InvalidParamsException(ErrorType.WORD_NOT_LESS_READ, w + "", readLength + "");
    }
    if (i > 0 && l > readLength - w) {
      throw new InvalidParamsException(ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, "-c", l + "", (readLength - w) + "");
    }
    if (i > readLength - w) {
      throw new InvalidParamsException(ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, "-b", i + "", (readLength - w) + "");
    }
    if (subs > readLength - w) {
      throw new InvalidParamsException(ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, "-a", subs + "", (readLength - w) + "");
    }
  }
}
