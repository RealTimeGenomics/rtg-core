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
import com.rtg.launcher.BuildCommon;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.HashingRegion;
import com.rtg.reader.FormatCli;
import com.rtg.reference.Sex;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;

/**
 * Flags common between different map commands
 */
public final class MapFlags {

  private static final int MAX_WORD_SIZE = 32; //Integer.valueOf(System.getProperty("rtg.max-word", "32"));
  private static final int MAX_INSERT_SIZE = 1000000000;

  private MapFlags() { }

  /** SDF input format */
  public static final String SDF_FORMAT = "sdf";
  /** FASTA input format */
  public static final String FASTA_FORMAT = "fasta";
  /** FASTQ input format */
  public static final String FASTQ_FORMAT = "fastq";
  /** SAM / BAM paired end input format */
  public static final String SAM_PE_FORMAT = "sam-pe";
  /** SAM / BAM single end input format */
  public static final String SAM_SE_FORMAT = "sam-se";
  static final String TSV_FORMAT = "tsv";
  static final String[] FORMAT_OPTIONS = {SDF_FORMAT, FASTA_FORMAT, FASTQ_FORMAT, SAM_SE_FORMAT, SAM_PE_FORMAT};
  /** Sanger quality format */
  public static final String SANGER_FORMAT = "sanger";
  /** Solexa quality format */
  public static final String SOLEXA_FORMAT = "solexa";
  /** Illumina quality format */
  public static final String ILLUMINA_FORMAT = "illumina";
  static final String[] QUALITY_FORMAT_OPTIONS = {SANGER_FORMAT, SOLEXA_FORMAT, ILLUMINA_FORMAT};

  static final String TEMPLATE_FLAG = "template";
  //private static final String INSERT_SIZE_FLAG = "insert-size";
  static final String THREAD_MULTIPLIER = "Xthread-multiplier";
  static final String X_LONG_READ = "Xlong-read";
  static final String OUTPUT_UNFILTERED = "all-hits";
  static final String OUTPUT_NULLFILTERED = "Xnull-filter";
  static final String READ_FREQUENCY_FLAG = "Xread-freq";
  static final String MIN_HITS_FLAG = "Xmin-hits";
  static final String NO_INMEMORY_TEMPLATE = CommonFlags.NO_INMEMORY_TEMPLATE;
  static final String OUTPUT_READ_NAMES_FLAG = "read-names";
  static final String LEGACY_CIGARS = "legacy-cigars";
  static final String BAM_FLAG = "bam";
  /** unify alignment outputs flag */
  public static final String DONT_UNIFY_FLAG = "no-merge";
  /** sam output flag */
  public static final String SAM_FLAG = "sam";
  /** Sex specification flag. */
  public static final String SEX_FLAG = "sex";
  static final String PEDIGREE_FLAG = "pedigree";
  static final String NO_CALIBRATION = "no-calibration";
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
    flags.getFlag(CommonFlags.INDELS_FLAG).setParameterDefault(1);
    flags.getFlag(CommonFlags.SUBSTITUTIONS_FLAG).setParameterDefault(1);

    flags.registerOptional('e', CommonFlags.MAX_ALIGNMENT_MISMATCHES, IntegerOrPercentage.class, CommonFlags.INT, "maximum mismatches for mappings in single-end mode (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf(NgsFilterParams.MAX_MATED_MISMATCH_THRESHOLD)).setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional('E', CommonFlags.UNMATED_MISMATCH_THRESHOLD, IntegerOrPercentage.class, CommonFlags.INT, "maximum mismatches for mappings of unmated results (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf(NgsFilterParams.MAX_UNMATED_MISMATCH_THRESHOLD)).setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional(CommonFlags.MATED_MISMATCH_THRESHOLD, IntegerOrPercentage.class, CommonFlags.INT, "maximum mismatches for mappings across mated results, alias for --" + CommonFlags.MAX_ALIGNMENT_MISMATCHES + " (as absolute value or percentage of read length)", IntegerOrPercentage.valueOf(NgsFilterParams.MAX_MATED_MISMATCH_THRESHOLD)).setCategory(CommonFlagCategories.REPORTING);

    flags.registerOptional(CommonFlags.TEMP_DIR, File.class, "DIR", "directory used for temporary files (Defaults to output directory)").setCategory(CommonFlagCategories.UTILITY);

    flags.registerOptional(LEGACY_CIGARS, "use legacy cigars in output").setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(OUTPUT_READ_NAMES_FLAG, "use read name in output instead of read id (Uses more RAM)").setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(ALIGNER_MODE_FLAG, AlignerMode.class, "STRING", "pick the aligner to use", AlignerMode.AUTO).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(SINGLE_INDEL_PENALTIES_FLAG, String.class, "STRING|FILE", "single indel penalty file", EditDistanceFactory.DEFAULT_SINGLE_INDEL_TABLE).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);

    //--X flags


    flags.registerOptional(THREAD_MULTIPLIER, Integer.class, CommonFlags.INT, "number of work chunks per thread", HashingRegion.DEFAULT_THREAD_MULTIPLIER).setCategory(CommonFlagCategories.UTILITY);
    //flags.registerOptional(INSERT_SIZE_FLAG, Integer.class, INT, "expected insert size for pairs");
    flags.registerOptional(CommonFlags.TEMP_FILES_COMPRESSED, Boolean.class, "BOOL", "gzip temporary SAM files", true).setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(NO_INMEMORY_TEMPLATE, "do not load the template in memory").setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(CommonFlags.FORCE_LONG_FLAG, "force the use of long read mode").setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional(SEX_FLAG, Sex.class, "sex", "sex of sample", null).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(MapFlags.PEDIGREE_FLAG, File.class, "file", "genome relationships pedigree containing sex of sample").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);

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

    flags.register(new Flag(null, MapFlags.ALIGNER_BAND_WIDTH_FACTOR_FLAG, "aligner indel band width scaling factor, fraction of read length allowed as an indel", 0, 1, Double.class, CommonFlags.FLOAT, 0.5, CommonFlagCategories.SENSITIVITY_TUNING));
  }

  /**
   * adds map related input and output flags
   * @param flags flags to set
   */
  public static void initInputOutputFlags(CFlags flags) {
    final Flag input = flags.registerOptional('i', CommonFlags.READS_FLAG, File.class, CommonFlags.SDF_OR_FILE, "input read set").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerRequired('o', CommonFlags.OUTPUT_FLAG, File.class, CommonFlags.DIR, BuildCommon.RESOURCE.getString("OUTPUT_DESC")).setCategory(CommonFlagCategories.INPUT_OUTPUT);
    CommonFlags.initTemplateFlag(flags);

    initInputFormatFlags(flags);
    flags.addRequiredSet(input);
  }

  /**
   * Initialise the input format flags for mapping commands
   * @param flags shared flags
   */
  public static void initInputFormatFlags(CFlags flags) {
    //mapx already has a -f flag
    final Flag formatFlag = flags.registerOptional('F', FormatCli.FORMAT_FLAG, String.class, "FORMAT", "input format for reads", SDF_FORMAT).setCategory(CommonFlagCategories.INPUT_OUTPUT);
    formatFlag.setParameterRange(FORMAT_OPTIONS);
    final Flag qualFormatFlag = flags.registerOptional('q', FormatCli.QUALITY_FLAG, String.class, "FORMAT", "format of quality data for fastq files (use sanger for Illumina 1.8+)").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    qualFormatFlag.setParameterRange(QUALITY_FORMAT_OPTIONS);
  }

  static void initPairedEndFormatFlags(CFlags flags)  {
    //IO format flags specific to paired end
    final Flag inputL = flags.registerOptional('l', FormatCli.LEFT_FILE_FLAG, File.class, "FILE", "left input file for FASTA/FASTQ paired end reads").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    final Flag inputR = flags.registerOptional('r', FormatCli.RIGHT_FILE_FLAG, File.class, "FILE", "right input file for FASTA/FASTQ paired end reads").setCategory(CommonFlagCategories.INPUT_OUTPUT);
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
    flags.registerOptional(DONT_UNIFY_FLAG, "output mated/unmated/unmapped alignments into separate SAM/BAM files").setCategory(CommonFlagCategories.INPUT_OUTPUT);
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
      final String format = flags.isSet(FormatCli.FORMAT_FLAG) ? flags.getValue(FormatCli.FORMAT_FLAG).toString().toLowerCase(Locale.getDefault()) : SDF_FORMAT;
      if (format.equals(SAM_SE_FORMAT) || format.equals(SAM_PE_FORMAT)) {
        flags.setParseMessage("Do not use left and right flags when using SAM single or paired end input.");
        return false;
      }
    }
    return true;
  }

  static boolean validateMapParams(CFlags flags) {
    final int a = (Integer) flags.getValue(CommonFlags.SUBSTITUTIONS_FLAG);
    final int b = (Integer) flags.getValue(CommonFlags.INDELS_FLAG);
    final int c = (Integer) flags.getValue(CommonFlags.INDEL_LENGTH_FLAG);
    if (flags.isSet(CommonFlags.WORDSIZE_FLAG)) {
      if (!CommonFlags.validateFlagBetweenValues(flags, CommonFlags.WORDSIZE_FLAG, 1, MAX_WORD_SIZE)) {
        return false;
      }
    }
    if (!CommonFlags.validateStepAndWordSize(flags)) {
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
    if (!CommonFlags.checkPercentRepeatFrequency(flags)) {
      return false;
    }

    if (!CommonFlags.validateThreads(flags)) {
      return false;
    }

    if (!CommonFlags.validateFlagBetweenValues(flags, CommonFlags.MAX_FRAGMENT_SIZE, 1, MAX_INSERT_SIZE)) {
      return false;
    }
    if (!CommonFlags.validateFlagBetweenValues(flags, CommonFlags.MIN_FRAGMENT_SIZE, 0, MAX_INSERT_SIZE)) {
      return false;
    }
    if (flags.isSet(CommonFlags.MATED_MISMATCH_THRESHOLD) && flags.isSet(CommonFlags.MAX_ALIGNMENT_MISMATCHES)) {
      flags.setParseMessage("--" + CommonFlags.MATED_MISMATCH_THRESHOLD + " is an alias for --" + CommonFlags.MAX_ALIGNMENT_MISMATCHES + ", please only set one");
      return false;
    }
    final IntegerOrPercentage maxMatedScore = (IntegerOrPercentage) flags.getValue(CommonFlags.MATED_MISMATCH_THRESHOLD);
    if (maxMatedScore.getValue(100) < 0) {
      Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + CommonFlags.MATED_MISMATCH_THRESHOLD, maxMatedScore + "", "0");
      return false;
    }
    final IntegerOrPercentage maxAlignmentScore = (IntegerOrPercentage) flags.getValue(CommonFlags.MAX_ALIGNMENT_MISMATCHES);
    if (maxAlignmentScore.getValue(100) < 0) {
      Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + CommonFlags.MAX_ALIGNMENT_MISMATCHES, maxAlignmentScore + "", "0");
      return false;
    }
    final IntegerOrPercentage maxUnMatedScore = (IntegerOrPercentage) flags.getValue(CommonFlags.UNMATED_MISMATCH_THRESHOLD);
    if (maxUnMatedScore.getValue(100) < 0) {
      Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + CommonFlags.UNMATED_MISMATCH_THRESHOLD, maxUnMatedScore + "", "0");
      return false;
    }

    if (!validatePenaltyFlags(flags)) {
      return false;
    }

    return true;
  }

  static boolean validateMapInputOutputParams(CFlags flags) {
    final String format = flags.isSet(FormatCli.FORMAT_FLAG) ? flags.getValue(FormatCli.FORMAT_FLAG).toString().toLowerCase(Locale.getDefault()) : SDF_FORMAT;
    final boolean sdf = format.equals(SDF_FORMAT);
    if (sdf) {
      if (!flags.isSet(CommonFlags.READS_FLAG)) {
        flags.setParseMessage("Must set --" + CommonFlags.READS_FLAG + " when --" + FormatCli.FORMAT_FLAG + " is set to " + SDF_FORMAT);
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
      flags.setParseMessage("--" + CommonFlags.MIN_FRAGMENT_SIZE + " requires --" + CommonFlags.READS_FLAG + " to be a paired end directory");
      return false;
    }
    if (flags.isSet(CommonFlags.MAX_FRAGMENT_SIZE) && !isPaired) {
      flags.setParseMessage("--" + CommonFlags.MAX_FRAGMENT_SIZE + " requires --" + CommonFlags.READS_FLAG + " to be a paired end directory");
      return false;
    }
    if (flags.isSet(CommonFlags.MATED_MISMATCH_THRESHOLD) && !isPaired) {
      flags.setParseMessage("--" + CommonFlags.MATED_MISMATCH_THRESHOLD + " requires --" + CommonFlags.READS_FLAG + " to be a paired end directory");
      return false;
    }
    if (flags.isSet(CommonFlags.UNMATED_MISMATCH_THRESHOLD) && !isPaired) {
      flags.setParseMessage("--" + CommonFlags.UNMATED_MISMATCH_THRESHOLD + " requires --" + CommonFlags.READS_FLAG + " to be a paired end directory");
      return false;
    }
    final int min = (Integer) flags.getValue(CommonFlags.MIN_FRAGMENT_SIZE);
    final int max = (Integer) flags.getValue(CommonFlags.MAX_FRAGMENT_SIZE);
    if (min > max) {
      flags.setParseMessage("--" + CommonFlags.MIN_FRAGMENT_SIZE + " must be less than --" + CommonFlags.MAX_FRAGMENT_SIZE);
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
    return CommonFlags.validateSexTemplateReference(flags, MapFlags.SEX_FLAG, MapFlags.PEDIGREE_FLAG, CommonFlags.TEMPLATE_FLAG);
  }

  static boolean validatePenaltyFlags(CFlags flags) {
    if (!CommonFlags.validateFlagBetweenValues(flags, MISMATCH_PENALTY_FLAG, 1, 100)) {
      return false;
    }
    if (!CommonFlags.validateFlagBetweenValues(flags, UNKNOWNS_PENALTY_FLAG, 0, 100)) {
      return false;
    }
    if (!CommonFlags.validateFlagBetweenValues(flags, GAP_OPEN_PENALTY_FLAG, 1, 100)) {
      return false;
    }
    if (!CommonFlags.validateFlagBetweenValues(flags, GAP_EXTEND_PENALTY_FLAG, 1, 100)) {
      return false;
    }
    if (!CommonFlags.validateFlagBetweenValues(flags, SOFT_CLIP_DISTANCE_FLAG, 0, 100)) {
      return false;
    }
    if (flags.isSet(ALIGNER_BAND_WIDTH_FACTOR_FLAG) && !CommonFlags.validateFlagBetweenValues(flags, ALIGNER_BAND_WIDTH_FACTOR_FLAG, 0.0f, 1.0f)) {
      return false;
    }
    if (!(AlignerMode.GENERAL == flags.getValue(ALIGNER_MODE_FLAG)) && (flags.isSet(MISMATCH_PENALTY_FLAG) || flags.isSet(UNKNOWNS_PENALTY_FLAG) || flags.isSet(GAP_OPEN_PENALTY_FLAG) || flags.isSet(GAP_EXTEND_PENALTY_FLAG))) {
      flags.setParseMessage("Penalty flags are only valid if --" + ALIGNER_MODE_FLAG + " is 'general'");
      return false;

    }

    return true;
  }
}
