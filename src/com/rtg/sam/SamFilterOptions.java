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
package com.rtg.sam;

import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;

import java.io.File;

import com.rtg.launcher.CommonFlags;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Flag;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;

/**
 * Constants and utility methods for command line flags for filtering SAM records.
 */
public final class SamFilterOptions {

  private SamFilterOptions() { }

  /** Flag name for filter on <code>IH</code> or <code>NH</code> attribute of SAM record. */
  public static final String MAX_HITS_FLAG = "max-hits";

  /** Used for a flag not having a single letter command. */
  public static final char NO_SINGLE_LETTER = '\0';

  private static final String HITS_DESC = "if set, ignore SAM records with an alignment count that exceeds this value";

  /**
   * Register flag for filtering on the number of hits.
   * @param flags flags to add into
   * @param singleLetter single letter code for option, or 0 for no single letter option
   * @return the flag
   */
  public static Flag registerMaxHitsFlag(final CFlags flags, final char singleLetter) {
    if (singleLetter != NO_SINGLE_LETTER) {
      return flags.registerOptional(singleLetter, MAX_HITS_FLAG, Integer.class, "int", HITS_DESC).setCategory(SENSITIVITY_TUNING);
    } else {
      return flags.registerOptional(MAX_HITS_FLAG, Integer.class, "int", HITS_DESC).setCategory(SENSITIVITY_TUNING);
    }
  }

  /** Flag name for filter on <code>MAPQ</code> field of SAM record. */
  public static final String MIN_MAPQ_FLAG = "min-mapq";

  private static final String MAPQ_DESC = "if set, ignore SAM records with MAPQ less than this value";

  /**
   * Register flag for filtering on the MAPQ
   *
   * @param flags flags to add into
   * @return the flag
   */
  public static Flag registerMinMapQFlag(final CFlags flags) {
    return flags.registerOptional(MIN_MAPQ_FLAG, Integer.class, "int", MAPQ_DESC).setCategory(SENSITIVITY_TUNING);
  }

  /** Flag name for filter of <code>AS</code> attribute of mated SAM records. */
  public static final String MAX_AS_MATED_FLAG = "max-as-mated";

  private static final String AS_MATED_DESC = "if set, ignore mated SAM records with an alignment score (AS attribute) that exceeds this value";

  /**
   * Register flag for mated <code>AS</code> filtering.
   * @param flags flags to add into
   * @param singleLetter single letter code for option, or 0 for no single letter option
   * @return the flag
   */
  public static Flag registerMaxASMatedFlag(final CFlags flags, final char singleLetter) {
    if (singleLetter != NO_SINGLE_LETTER) {
      return flags.registerOptional(singleLetter, MAX_AS_MATED_FLAG, IntegerOrPercentage.class, CommonFlags.INT, AS_MATED_DESC).setCategory(SENSITIVITY_TUNING);
    } else {
      return flags.registerOptional(MAX_AS_MATED_FLAG, IntegerOrPercentage.class, CommonFlags.INT, AS_MATED_DESC).setCategory(SENSITIVITY_TUNING);
    }
  }

  /** Flag name for filter of <code>AS</code> attribute of unmated SAM records. */
  public static final String MAX_AS_UNMATED_FLAG = "max-as-unmated";

  private static final String AS_UNMATED_DESC = "if set, ignore unmated SAM records with an alignment score (AS attribute) that exceeds this value";

  /**
   * Register flag for unmated <code>AS</code> filtering.
   * @param flags flags to add into
   * @param singleLetter single letter code for option, or 0 for no single letter option
   * @return the flag
   */
  public static Flag registerMaxASUnmatedFlag(final CFlags flags, final char singleLetter) {
    if (singleLetter != NO_SINGLE_LETTER) {
      return flags.registerOptional(singleLetter, MAX_AS_UNMATED_FLAG, IntegerOrPercentage.class, CommonFlags.INT, AS_UNMATED_DESC).setCategory(SENSITIVITY_TUNING);
    } else {
      return flags.registerOptional(MAX_AS_UNMATED_FLAG, IntegerOrPercentage.class, CommonFlags.INT, AS_UNMATED_DESC).setCategory(SENSITIVITY_TUNING);
    }
  }

  /** Flag name for filtering out mated results. */
  public static final String EXCLUDE_MATED_FLAG = "exclude-mated";

  private static final String EXCLUDE_MATED_DESC = "exclude all mated SAM records";

  /**
   * Register flag for excluding mated results.
   *
   * @param flags flags to add into
   * @return the flag
   */
  public static Flag registerExcludeMatedFlag(final CFlags flags) {
    return flags.registerOptional(EXCLUDE_MATED_FLAG, EXCLUDE_MATED_DESC).setCategory(SENSITIVITY_TUNING);
  }

  /** Flag name for filtering out unmated results. */
  public static final String EXCLUDE_UNMATED_FLAG = "exclude-unmated";

  private static final String EXCLUDE_UNMATED_DESC = "exclude all unmated SAM records";

  /**
   * Register flag for excluding unmated results.
   *
   * @param flags flags to add into
   * @return the flag
   */
  public static Flag registerExcludeUnmatedFlag(final CFlags flags) {
    return flags.registerOptional(EXCLUDE_UNMATED_FLAG, EXCLUDE_UNMATED_DESC).setCategory(SENSITIVITY_TUNING);
  }

  /** Flag name to enable filtering out unmapped results. */
  public static final String EXCLUDE_UNMAPPED_FLAG = "exclude-unmapped";

  private static final String EXCLUDE_UNMAPPED_DESC = "exclude all unmapped SAM records";

  /**
   * Register flag for excluding unmapped results.
   *
   * @param flags flags to add into
   * @return the flag
   */
  public static Flag registerExcludeUnmappedFlag(final CFlags flags) {
    return flags.registerOptional(EXCLUDE_UNMAPPED_FLAG, EXCLUDE_UNMAPPED_DESC).setCategory(SENSITIVITY_TUNING);
  }

  /** Flag name to enable filtering out unplaced results. */
  public static final String EXCLUDE_UNPLACED_FLAG = "exclude-unplaced";

  private static final String EXCLUDE_UNPLACED_DESC = "exclude all SAM records with no alignment position";

  /**
   * Register flag for excluding unmapped results.
   *
   * @param flags flags to add into
   * @return the flag
   */
  public static Flag registerExcludeUnplacedFlag(final CFlags flags) {
    return flags.registerOptional(EXCLUDE_UNPLACED_FLAG, EXCLUDE_UNPLACED_DESC).setCategory(SENSITIVITY_TUNING);
  }

  /** Flag name to enable filtering out duplicate results. */
  public static final String EXCLUDE_DUPLICATES_FLAG = "exclude-duplicates";

  private static final String EXCLUDE_DUPLICATES_DESC = "exclude all SAM records flagged as a PCR or optical duplicate";

  /**
   * Register flag for excluding duplicate results.
   *
   * @param flags flags to add into
   * @return the flag
   */
  public static Flag registerExcludeDuplicatesFlag(final CFlags flags) {
    return flags.registerOptional(EXCLUDE_DUPLICATES_FLAG, EXCLUDE_DUPLICATES_DESC).setCategory(SENSITIVITY_TUNING);
  }

  // this options differs from the above as it refers to the mechanism of detecting and removing duplicates on the fly, as well as looking at the SAM records FLAG field
  /** keep duplicates flag constant */
  public static final String KEEP_DUPLICATES_FLAG = "keep-duplicates";

  private static final String KEEP_DUPLICATES_DESC = "don't detect and filter duplicate reads based on mapping position";

  /**
   * Register flag for keeping duplicate results.
   * @param flags flags to add into
   * @return the flag
   */
  public static Flag registerKeepDuplicatesFlag(final CFlags flags) {
    return flags.registerOptional(KEEP_DUPLICATES_FLAG, KEEP_DUPLICATES_DESC).setCategory(SENSITIVITY_TUNING);
  }

  /** Flag name for filtering out uninteresting sequences. */
  public static final String RESTRICTION_FLAG = "region";

  private static final String RESTRICTION_DESC = "if set, only process SAM records within the specified range. The format is one of <template_name>, <template_name>:start-end or <template_name>:start+length";

  /**
   * Register flag for restricting records to be processed.
   *
   * @param flags flags to add into
   * @return the flag
   */
  public static Flag registerRestrictionFlag(final CFlags flags) {
    return flags.registerOptional(RESTRICTION_FLAG, String.class, "string", RESTRICTION_DESC).setCategory(SENSITIVITY_TUNING);
  }

  /** Flag for a bed region restriction */
  public static final String BED_REGIONS_FLAG = "bed-regions";

  /**
   * Register flag for restricting records to be processed.
   *
   * @param flags flags to add into
   * @return the flag
   */
  public static Flag registerBedRestrictionFlag(final CFlags flags) {
    return flags.registerOptional(BED_REGIONS_FLAG, File.class, "FILE", "BED file containing regions to process").setCategory(SENSITIVITY_TUNING);
  }


  private static final String REQUIRE_FLAGS = "require-flags";
  private static final String FILTER_FLAGS = "filter-flags";

  /**
   * Register flags for restricting directly based on flags.
   * @param flags the flags to add in to
   */
  public static void registerMaskFlags(CFlags flags) {
    flags.registerOptional('f', REQUIRE_FLAGS, Integer.class, "INT", "decimal mask indicating SAM FLAG bits that must be set for the record").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('F', FILTER_FLAGS, Integer.class, "INT", "decimal mask indicating SAM FLAG bits that must not be set for the record").setCategory(SENSITIVITY_TUNING);
  }

  /**
   * Validate the filter flag input, will produce the appropriate
   * diagnostic error when validation fails.
   *
   * @param flags the flags object to check
   * @param allowUnmappedOnly true if the user is allowed to exclude both mated and unmated alignments
   * @return true if the provided flags are valid, false otherwise
   */
  public static boolean validateFilterFlags(final CFlags flags, boolean allowUnmappedOnly) {
    if (flags.isSet(MAX_HITS_FLAG)) {
      final int maxHits = (Integer) flags.getValue(MAX_HITS_FLAG);
      if (maxHits < 1) {
        Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + MAX_HITS_FLAG, "" + maxHits, "1");
        return false;
      }
    }
    if (flags.isSet(MIN_MAPQ_FLAG)) {
      final int minMapQ = (Integer) flags.getValue(MIN_MAPQ_FLAG);
      if (minMapQ < 1) {
        Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + MIN_MAPQ_FLAG, "" + minMapQ, "1");
        return false;
      }
    }
    if (flags.isSet(MAX_AS_MATED_FLAG)) {
      final IntegerOrPercentage maxMated = (IntegerOrPercentage) flags.getValue(MAX_AS_MATED_FLAG);
      if (maxMated.getValue(100) < 0) {
        Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + MAX_AS_MATED_FLAG, "" + maxMated, "0");
        return false;
      }
    }
    if (flags.isSet(MAX_AS_UNMATED_FLAG)) {
      final IntegerOrPercentage maxUnmated = (IntegerOrPercentage) flags.getValue(MAX_AS_UNMATED_FLAG);
      if (maxUnmated.getValue(100) < 0) {
        Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + MAX_AS_UNMATED_FLAG, "" + maxUnmated, "0");
        return false;
      }
    }
    if (!(allowUnmappedOnly || flags.checkNand(EXCLUDE_MATED_FLAG, EXCLUDE_UNMATED_FLAG))) {
      return false;
    }

    if (flags.isSet(FILTER_FLAGS) && flags.isSet(REQUIRE_FLAGS)) {
      final int unset = (Integer) flags.getValue(FILTER_FLAGS);
      final int set = (Integer) flags.getValue(REQUIRE_FLAGS);
      final int badFlags = unset & set;
      if (badFlags != 0) {
        flags.setParseMessage("--" + FILTER_FLAGS + " and --" + REQUIRE_FLAGS + " have conflicting values. Flags in common: " + badFlags);
      }
    }

    if (flags.isSet(RESTRICTION_FLAG)) {
      final String region = (String) flags.getValue(RESTRICTION_FLAG);
      if (!RegionRestriction.validateRegion(region)) {
        flags.setParseMessage("The value \"" + region + "\" for \"--" + RESTRICTION_FLAG + "\" is malformed.");
        return false;
      }
    }
    if (flags.isSet(BED_REGIONS_FLAG)) {
      final File bedRegionsFile = (File) flags.getFlag(BED_REGIONS_FLAG).getValue();
      if (!bedRegionsFile.exists()) {
        Diagnostic.error(ErrorType.FILE_NOT_FOUND,
          "The specified file, \"" + bedRegionsFile.getPath() + "\", does not exist.");
        return false;
      }
    }
    if (flags.isSet(BED_REGIONS_FLAG) && flags.isSet(SamFilterOptions.RESTRICTION_FLAG)) {
      flags.setParseMessage("Can not specify both --" + BED_REGIONS_FLAG + " and --" + SamFilterOptions.RESTRICTION_FLAG);
      return false;
    }

    if (flags.getFlag(KEEP_DUPLICATES_DESC) != null
        && flags.getFlag(EXCLUDE_DUPLICATES_FLAG) != null) {
      throw new RuntimeException("Cannot have registered flags for both include and exclude duplicates");
    }
    return true;
  }

  /**
   * Build parameters from flags.
   * @param flags command line flags
   * @return parameters
   */
  public static SamFilterParams.SamFilterParamsBuilder makeFilterParamsBuilder(final CFlags flags) {
    final SamFilterParams.SamFilterParamsBuilder builder = SamFilterParams.builder();
    if (flags.isSet(MAX_HITS_FLAG)) {
      builder.maxAlignmentCount((Integer) flags.getValue(MAX_HITS_FLAG));
    }
    if (flags.isSet(MIN_MAPQ_FLAG)) {
      builder.minMapQ((Integer) flags.getValue(MIN_MAPQ_FLAG));
    }
    if (flags.isSet(MAX_AS_MATED_FLAG)) {
      final IntegerOrPercentage matedAS = (IntegerOrPercentage) flags.getValue(MAX_AS_MATED_FLAG);
      builder.maxMatedAlignmentScore(matedAS);
      if (flags.isSet(MAX_AS_UNMATED_FLAG)) {
        final IntegerOrPercentage unmatedAS = (IntegerOrPercentage) flags.getValue(MAX_AS_UNMATED_FLAG);
        if (unmatedAS.compareTo(matedAS) > 0) {
          Diagnostic.warning("--" + MAX_AS_UNMATED_FLAG + " should not be greater than --" + MAX_AS_MATED_FLAG);
        }
      }
    }
    if (flags.isSet(MAX_AS_UNMATED_FLAG)) {
      builder.maxUnmatedAlignmentScore((IntegerOrPercentage) flags.getValue(MAX_AS_UNMATED_FLAG));
    }
    builder.excludeMated(flags.isSet(EXCLUDE_MATED_FLAG));
    builder.excludeUnmated(flags.isSet(EXCLUDE_UNMATED_FLAG));
    builder.excludeUnmapped(flags.isSet(SamFilterOptions.EXCLUDE_UNMAPPED_FLAG));
    builder.excludeUnplaced(flags.isSet(SamFilterOptions.EXCLUDE_UNPLACED_FLAG));
    if (flags.isSet(FILTER_FLAGS)) {
      builder.requireUnsetFlags((Integer) flags.getValue(FILTER_FLAGS));
    }
    if (flags.isSet(REQUIRE_FLAGS)) {
      builder.requireSetFlags((Integer) flags.getValue(REQUIRE_FLAGS));
    }

    // Some tools want inclusion by default and some want exclusion by default
    if (flags.getFlag(KEEP_DUPLICATES_FLAG) != null) {
      builder.findAndRemoveDuplicates(!flags.isSet(KEEP_DUPLICATES_FLAG));
      builder.excludeDuplicates(!flags.isSet(KEEP_DUPLICATES_FLAG));
    } else if (flags.getFlag(EXCLUDE_DUPLICATES_FLAG) != null) {
      builder.excludeDuplicates(flags.isSet(EXCLUDE_DUPLICATES_FLAG));
    }

    if (flags.isSet(RESTRICTION_FLAG)) {
      builder.restriction((String) flags.getValue(RESTRICTION_FLAG));
    }
    if (flags.isSet(BED_REGIONS_FLAG)) {
      builder.bedRegionsFile((File) flags.getValue(BED_REGIONS_FLAG));
    }

    return builder;
  }
}
