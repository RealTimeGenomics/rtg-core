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
package com.rtg.launcher;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Flag;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.AdjustableGZIPOutputStream;

/**
 * A global store of command line flags to use for highly experimental options
 * <p>
 *   NOTE care should be taken to ensure flags have been set by the time your code accessing the values is initialised
 *   (i.e. they should not get values during the static initialisation of a <code>CLI</code> class, since this could
 *   cause them to access the values before the <code>CFlags.setFlags(args)</code> call is made. However other things
 *   can also cause a classes static initialization to start)
 * </p>
 */
public final class GlobalFlags {

  private GlobalFlags() { }
  //NOTE flag definitions should be placed here, to prevent circular dependencies
  //NOTE read class note

  //Map
  /** keep temporary files from mapping run instead of deleting them */
  public static final String MAP_KEEP_TEMPORARY_FILES = "com.rtg.map.keep-temporary-files";

  //Species
  /** Threshold for termination in Species - when smoothed estimate of remaining changes to L drops below this then terminate solving loop.  */
  public static final String SPECIES_LTERMINATION_FLAG = "com.rtg.species.ltermination";
  /** Test early termination for p-values */
  public static final String SPECIES_TERMINATION_TARGET_FLAG = "com.rtg.species.ltermination-target";

  //Edit distance factory
  /** Specify how many reads to log */
  public static final String EDIT_DIST_LOGGING_AMOUNT_FLAG = "com.rtg.alignment.EditDistanceFactory.logging-amount";
  /** Enable the heuristic aligners (faster, but some lower quality alignments are produced) */
  public static final String EDIT_DIST_HEURISTIC_ALIGNERS_FLAG = "com.rtg.alignment.EditDistanceFactory.heuristic-aligners";
  /** Only use the Gotoh aligner (disable all others) */
  public static final String EDIT_DIST_GOTOH_ONLY_FLAG = "com.rtg.alignment.EditDistanceFactory.gotoh-only";
  /** Only use the SingleIndelEditDistance aligner (disable all others) */
  public static final String EDIT_DIST_SINGLE_INDEL_ONLY_FLAG = "com.rtg.alignment.EditDistanceFactory.single-indel-only";
  /** Only use the SingleIndelSeededEditDistance aligner (disable all others) */
  public static final String EDIT_DIST_SINGLE_INDEL_SEEDED_ONLY_FLAG = "com.rtg.alignment.EditDistanceFactory.single-indel-seeded-only";

  /** If set, load insertion / deletion penalty distribution from the supplied file */
  public static final String EDIT_DIST_INDEL_TABLE_FLAG = "com.rtg.alignment.SingleIndelEditDistance.penalties-file";

  /** True to log alignment score histogram */
  public static final String EDIT_DIST_LOG_AS_HISTOGRAM_FLAG = "com.rtg.alignment.UnidirectionalPrioritisedEditDistance.log-as-histo";

  /** Dump alignment stats upon closing a temp file writer */
  public static final String TEMP_FILES_DUMP_ALIGN_STATS_FLAG = "com.rtg.ngs.tempstage.AbstractTempFileWriter.dump-alignment-stats";

  /** Allow SAM file loading when sam header is not coordinate sorted */
  public static final String SAM_IGNORE_SORT_ORDER_FLAG = "com.rtg.sam.ignore-header-sortorder";

  /** Allow fallback to a slower alternative when reading non-indexed SAM files with region restrictions requested */
  public static final String SAM_ALLOW_FALLBACK_FOR_NON_INDEXED_REGIONS = "com.rtg.sam.allow-region-fallback";

  /** Maximum number of hits at a given position in the sliding window collector */
  //see bug #1476 for consequences of this on larger datasets
  public static final String SLIDING_WINDOW_MAX_HITS_PER_POS_FLAG = "com.rtg.pairedend.SlidingWindow.max-hits-per-position";

  /** Default total length of all inserts/deletes allowed in reasonably short reads */
  public static final String DEFAULT_INDEL_LENGTH_FLAG = "com.rtg.util.default-indel-length";

  /** If true, the population command will fall back to using forward backward when disagreeing calls are encountered (currently slow for large pops) */
  public static final String FAMILY_CALLER_FALLBACK_FLAG = "com.rtg.variant.bayes.multisample.FamilyCaller.fb-fallback";

  /** The maximum number of hypotheses that can comfortably be handled by the complex caller */
  public static final String COMPLEX_CALLER_MAX_HYPOTH_FLAG = "com.rtg.variant.bayes.multisample.ComplexCaller.max-hypoth";

  /** Do calling aware of ploidy in PAR regions */
  public static final String CALLER_PAR_AWARE = "com.rtg.variant.par-aware";

  /** Complex region extraction include indel lengths in interesting separation */
  public static final String COMPLEX_REGION_INDEL_EXTENSION = "com.rtg.variant.region-indel-extension";

  /** Complex region extraction maximum unit size looked for by SimpleRepeatMeasurer, e.g. 3-mer repeats */
  public static final String COMPLEX_REGION_SIMPLE_REPEAT_LIMIT = "com.rtg.variant.region-simple-repeat-limit";

  /** Complex region extraction simple repeat implementation */
  public static final String COMPLEX_REGION_SIMPLE_REPEAT_IMPL = "com.rtg.variant.region-simple-repeat-impl";

  /** Variant caller min depth for call-at-N triggering */
  public static final String CALLER_N_MIN_DEPTH = "com.rtg.variant.n-min-depth";

  /** Prior for random hypothesis when using random hypotheses. */
  public static final String ENTROPY_RANDOM_PRIOR_FLAG = "com.rtg.variant.bayes.EntropyRandom.random-prior";

  //Assembler
  /** If more than this many hits are seen at a position, skip them all. */
  public static final String ASSEMBLER_MAX_HITS_PER_START_POS_FLAG = "com.rtg.assembler.maxhits";
  /** Number of deviations to apply to insert distributions. */
  public static final String ASSEMBLER_INSERT_DEVIATIONS_FLAG = "com.rtg.assembler.insertdeviations";

  /** Should the random tree builder push missing values down during the build process. */
  public static final String TRAIN_ON_MISSING_VALUES = "com.rtg.ml.train-on-missing";

  /** Specify the maximum number of simultaneous paths before vcfeval skips a region */
  public static final String VCFEVAL_MAX_PATHS = "com.rtg.variant.eval.max-paths";

  /** Turn on alternate ROC slope calculation */
  public static final String ALTERNATE_ROC_SLOPE_CALCULATION = "com.rtg.variant.eval.RocSlope.alt-roc-slope";

  /** Level of BAM compression to use during recalibration (probably also works for SAM merge). */
  public static final String GZIP_LEVEL = "com.rtg.calibrate.Recalibrate.gzip-level";


  private static CFlags sFlags;
  private static final Set<String> ACCESSED_FLAGS = new HashSet<>();
  /** Category for global flags */
  public static final String CATEGORY = "Highly Experimental";
  private static final ArrayList<Flag> FLAGS = new ArrayList<>();

  static {

    // Metagenomics
    registerFlag(SPECIES_LTERMINATION_FLAG, Double.class, 0.1);
    registerFlag(SPECIES_TERMINATION_TARGET_FLAG, Double.class, 0.01);

    //Edit distance factory
    registerFlag(EDIT_DIST_LOGGING_AMOUNT_FLAG, Integer.class, 0);
    registerFlag(EDIT_DIST_HEURISTIC_ALIGNERS_FLAG, Boolean.class, true);
    registerFlag(EDIT_DIST_GOTOH_ONLY_FLAG);
    registerFlag(EDIT_DIST_SINGLE_INDEL_ONLY_FLAG);
    registerFlag(EDIT_DIST_SINGLE_INDEL_SEEDED_ONLY_FLAG);
    registerFlag(EDIT_DIST_INDEL_TABLE_FLAG, String.class, "");
    registerFlag(EDIT_DIST_LOG_AS_HISTOGRAM_FLAG);

    registerFlag(TEMP_FILES_DUMP_ALIGN_STATS_FLAG);
    registerFlag(MAP_KEEP_TEMPORARY_FILES);
    registerFlag(SLIDING_WINDOW_MAX_HITS_PER_POS_FLAG, Integer.class, 1000);
    registerFlag(SAM_IGNORE_SORT_ORDER_FLAG);
    registerFlag(SAM_ALLOW_FALLBACK_FOR_NON_INDEXED_REGIONS);

    // Aligners / all-paths
    registerFlag(DEFAULT_INDEL_LENGTH_FLAG, Integer.class, 7);

    registerFlag(ASSEMBLER_MAX_HITS_PER_START_POS_FLAG, Integer.class, 5);
    registerFlag(ASSEMBLER_INSERT_DEVIATIONS_FLAG, Integer.class, 4);

    // Complex caller
    registerFlag(COMPLEX_CALLER_MAX_HYPOTH_FLAG, Integer.class, 20);
    registerFlag(CALLER_PAR_AWARE, Boolean.class, false);
    registerFlag(COMPLEX_REGION_INDEL_EXTENSION, Boolean.class, false);
    registerFlag(COMPLEX_REGION_SIMPLE_REPEAT_LIMIT, Integer.class, 3);
    registerFlag(COMPLEX_REGION_SIMPLE_REPEAT_IMPL, String.class, "default");

    // Misc calling
    registerFlag(CALLER_N_MIN_DEPTH, Integer.class, 5);
    registerFlag(FAMILY_CALLER_FALLBACK_FLAG, Boolean.class, false);
    registerFlag(ENTROPY_RANDOM_PRIOR_FLAG, Double.class, 0.0);

    // AVR, training on missing instances increases time and experience indicates is a bad idea
    // when there are lots of missing values.
    registerFlag(TRAIN_ON_MISSING_VALUES, Boolean.class, false);

    registerFlag(VCFEVAL_MAX_PATHS, Integer.class, 5000);
    registerFlag(ALTERNATE_ROC_SLOPE_CALCULATION);
    registerFlag(GZIP_LEVEL, Integer.class, AdjustableGZIPOutputStream.DEFAULT_GZIP_LEVEL);
  }

  private static final CFlags DEFAULT_FLAGS = new CFlags();
  static { //this ensures default values are available in tests and classes which don't use <code>CLI</code> framework
    registerExperimentalFlags(DEFAULT_FLAGS);
  }

  static void registerFlag(String name, Class<?> type, Object def) {
    if (type != null && def == null) {
      throw new IllegalArgumentException("Default value must be non-null for experimental flags with a type");
    }
    FLAGS.add(new Flag(null, "XX" + name, "N/A", 0, 1, type, "N/A", def, CATEGORY));
  }

  private static void registerFlag(String name) {
    registerFlag(name, null, null);
  }

  static void registerExperimentalFlags(CFlags flags) {
    final String[] cat = flags.getCategories();
    if (cat != null) {
      final String[] copy = Arrays.copyOf(cat, cat.length + 1);
      copy[copy.length - 1] = CATEGORY;
      flags.setCategories(flags.getHelpCategory(), copy);
    }
    for (final Flag flag : FLAGS) {
      flags.register(flag);
    }
    sFlags = flags;
    resetAccessedStatus();
  }

  /**
   * @param name name of flag (without <code>XX</code> prefix)
   * @return flag object
   */
  public static Flag getFlag(String name) {
    final String innerName = "XX" + name;
    ACCESSED_FLAGS.add(innerName);
    return sFlags.getFlag(innerName);
  }

  /**
   * @param name name of flag (without <code>XX</code> prefix)
   * @return true if flag is set
   */
  public static boolean isSet(String name) {
    return getFlag(name).isSet();
  }

  /**
   * @param name name of flag (without <code>XX</code> prefix)
   * @return string representation of flag value
   */
  public static String getStringValue(String name) {
    return (String) getFlag(name).getValue();
  }

  /**
   * @param name name of flag (without <code>XX</code> prefix)
   * @return boolean representation of flag value
   */
  public static boolean getBooleanValue(String name) {
    return (Boolean) getFlag(name).getValue();
  }

  /**
   * @param flagName name of flag (without <code>xx</code> prefix)
   * @return int representation of flag value
   */
  public static int getIntegerValue(String flagName) {
    return (Integer) getFlag(flagName).getValue();
  }

  /**
   * checks flags haven't been accessed yet, prints a warning if they have
   * @return true if no flags have been accessed
   */
  public static boolean initialAccessCheck() {
    boolean bad = false;
    for (final String flag : ACCESSED_FLAGS) {
      Diagnostic.warning("Flag: " + flag + " is accessed before flag registration");
      bad = true;
    }
    resetAccessedStatus();
    return !bad;
  }

  /**
   * @return true iff all set flags were also accessed since registration
   */
  public static boolean finalAccessCheck() {
    for (final Flag f : FLAGS) {
      if (f.isSet() && !ACCESSED_FLAGS.contains(f.getName())) {
        Diagnostic.warning("Flag: " + f.getName() + " is set but never accessed.");
        resetAccessedStatus();
        return false;
      }
    }
    resetAccessedStatus();
    return true;
  }

  /**
   * Unsets list of flags which have been accessed. Use this if your test is failing due to another test using a global flag
   */
  public static void resetAccessedStatus() {
    ACCESSED_FLAGS.clear();
  }
}
