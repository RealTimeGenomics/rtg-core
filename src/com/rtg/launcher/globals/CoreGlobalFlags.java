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
package com.rtg.launcher.globals;

import java.util.List;

import com.reeltwo.jumble.annotations.JumbleIgnore;
import com.rtg.util.cli.Flag;
import com.rtg.util.io.AdjustableGZIPOutputStream;

/**
 * Experimental flags for public release
 */
@JumbleIgnore
public class CoreGlobalFlags extends GlobalFlagsInitializer {

  /** Threshold for termination in Species - when smoothed estimate of remaining changes to L drops below this then terminate solving loop.  */
  public static final String SPECIES_LTERMINATION_FLAG = "com.rtg.species.ltermination";
  /** Test early termination for p-values */
  public static final String SPECIES_TERMINATION_TARGET_FLAG = "com.rtg.species.ltermination-target";
  /** Allow shut off of principle component analysis in similarity. */
  public static final String SIMILARITY_PCA_FLAG = "com.rtg.similarity.pca";
  /** Specify how many reads to log */
  public static final String EDIT_DIST_LOGGING_AMOUNT_FLAG = "com.rtg.alignment.EditDistanceFactory.logging-amount";
  /** Enable the heuristic aligners (faster, but some lower quality alignments are produced) */
  public static final String EDIT_DIST_HEURISTIC_ALIGNERS_FLAG = "com.rtg.alignment.EditDistanceFactory.heuristic-aligners";
  /** Only use the Gotoh aligner (disable all others) */
  public static final String EDIT_DIST_GOTOH_ONLY_FLAG = "com.rtg.alignment.EditDistanceFactory.gotoh-only";
  /** Only use the <code>SingleIndelSeededEditDistance</code> aligner (disable all others) */
  public static final String EDIT_DIST_SINGLE_INDEL_SEEDED_ONLY_FLAG = "com.rtg.alignment.EditDistanceFactory.single-indel-seeded-only";
  /** If set, load insertion / deletion penalty distribution from the supplied file */
  public static final String EDIT_DIST_INDEL_TABLE_FLAG = "com.rtg.alignment.SingleIndelEditDistance.penalties-file";
  /** True to log alignment score histogram */
  public static final String EDIT_DIST_LOG_AS_HISTOGRAM_FLAG = "com.rtg.alignment.UnidirectionalPrioritisedEditDistance.log-as-histo";
  /** Dump alignment stats upon closing a temp file writer */
  public static final String TEMP_FILES_DUMP_ALIGN_STATS_FLAG = "com.rtg.ngs.tempstage.AbstractTempFileWriter.dump-alignment-stats";
  /** keep temporary files from mapping run instead of deleting them */
  public static final String MAP_KEEP_TEMPORARY_FILES = "com.rtg.map.keep-temporary-files";
  /** Maximum number of hits at a given position in the sliding window collector */
  //see bug #1476 for consequences of this on larger datasets
  public static final String SLIDING_WINDOW_MAX_HITS_PER_POS_FLAG = "com.rtg.pairedend.SlidingWindow.max-hits-per-position";
  /** Default total length of all inserts/deletes allowed in reasonably short reads */
  public static final String DEFAULT_INDEL_LENGTH_FLAG = "com.rtg.util.default-indel-length";
  /** If more than this many hits are seen at a position, skip them all. */
  public static final String ASSEMBLER_MAX_HITS_PER_START_POS_FLAG = "com.rtg.assembler.maxhits";
  /** Number of deviations to apply to insert distributions. */
  public static final String ASSEMBLER_INSERT_DEVIATIONS_FLAG = "com.rtg.assembler.insertdeviations";
  /** The maximum number of hypotheses that can comfortably be handled by the complex caller */
  public static final String COMPLEX_CALLER_MAX_HYPOTH_FLAG = "com.rtg.variant.bayes.multisample.ComplexCaller.max-hypoth";
  /** Complex region extraction include indel lengths in interesting separation */
  public static final String COMPLEX_REGION_INDEL_EXTENSION = "com.rtg.variant.region-indel-extension";
  /** Complex region extraction maximum unit size looked for by <code>SimpleRepeatMeasurer</code>, e.g. 3-mer repeats */
  public static final String COMPLEX_REGION_SIMPLE_REPEAT_LIMIT = "com.rtg.variant.region-simple-repeat-limit";
  /** Complex region extraction simple repeat implementation */
  public static final String COMPLEX_REGION_SIMPLE_REPEAT_IMPL = "com.rtg.variant.region-simple-repeat-impl";
  /** Dump the non SNP signals that are used for triggering complex calls */
  public static final String DUMP_COMPLEX_TRIGGER_SIGNALS = "com.rtg.variant.dump-complex-trigger-signals";
  /** Use soft clips to trigger complex calling */
  public static final String SOFT_CLIP_COMPLEX_TRIGGER = "com.rtg.variant.soft-clip-complex-trigger";
  /** Print the details of complex evidence into screen. */
  public static final String COMPLEX_EVIDENCE_DETAILS = "com.rtg.variant.complex.complex-evidence-details";
  /** Print the details of complex hypotheses into screen. */
  public static final String COMPLEX_HYPOTHESIS_DETAILS = "com.rtg.variant.complex.complex-hypothesis-details";
  /** Use new method for computing complex hypothesis priors */
  public static final String COMPLEX_HYPOTHESIS_NEW_PRIORS = "com.rtg.variant.complex.complex-hypothesis-new-priors";
  /** Adjust complex diploid priors using genome priors and other empirical tweaks */
  public static final String COMPLEX_HYPOTHESIS_ADJUST_PRIORS = "com.rtg.variant.complex.complex-hypothesis-adjust-priors";
  /** If true, all-paths should attempt to use unrolled CG read, otherwise use the flattened representation */
  public static final String COMPLEX_CALLER_UNROLL_CG_FLAG = "com.rtg.variant.bayes.EvidenceComplex.unroll-cg";
  /** Variant caller min depth for call-at-N triggering */
  public static final String CALLER_N_MIN_DEPTH = "com.rtg.variant.n-min-depth";
  /** If true, the population command will fall back to using forward backward when disagreeing calls are encountered (currently slow for large pops) */
  public static final String FAMILY_CALLER_FALLBACK_FLAG = "com.rtg.variant.bayes.multisample.FamilyCaller.fb-fallback";
  /** If true perform early exit of family caller when enough precision is reached */
  public static final String FAMILY_CALLER_SORTED_HYPOTHESES = "com.rtg.variant.bayes.multisample.FamilyPosterior.sorted-hypotheses";
  /** Treat bases with phred below the minimum base quality as quality 2*/
  public static final String MIN_BASE_QUALITY_AS_TWO = "com.rtg.variant.VariantAlignmentRecord.min-bq-2";
  /** Minimum count required in a calibration covariate set before it will be used */
  public static final String QUALITY_CALIBRATION_MIN_EVIDENCE = "com.rtg.variant.quality-calibration-min-evidence";
  /** Use covariate intersection calibration method */
  public static final String QUALITY_CALIBRATION_COVARIATE_INTERSECTION = "com.rtg.variant.quality-calibration-covariate-intersection";
  /** Use Dirichlet rather than multinomial allele balance calculation. */
  public static final String TUMOR_ALLELE_BALANCE = "com.rtg.variant.bayes.multisample.cancer.allele-balance-type";
  /** Should the random tree builder push missing values down during the build process. */
  public static final String AVR_TRAIN_ON_MISSING_VALUES = "com.rtg.ml.train-on-missing";
  /** Allow prediction to continue even if the VCF does not declare all the attributes of the model. */
  public static final String AVR_ALLOW_UNDECLARED_ATTRIBUTES = "com.rtg.avr.allow-undeclared";
  /** Level of BAM compression to use during recalibration (probably also works for SAM merge). */
  public static final String GZIP_LEVEL = "com.rtg.calibrate.Recalibrate.gzip-level";
  /** Allow fallback to a slower alternative when reading non-indexed SAM files with region restrictions requested */
  public static final String SAM_ALLOW_FALLBACK_FOR_NON_INDEXED_REGIONS = "com.rtg.sam.allow-region-fallback";
  /** Use code page 437 characters for deletes*/
  public static final String CP437_DELETES = "com.rtg.visualization.cp437-deletes";

  CoreGlobalFlags(List<Flag> flags) {
    super(flags);
  }

  @Override
  public void registerFlags() {
    // Metagenomics
    registerFlag(SPECIES_LTERMINATION_FLAG, Double.class, 0.1);
    registerFlag(SPECIES_TERMINATION_TARGET_FLAG, Double.class, 0.01);
    registerFlag(SIMILARITY_PCA_FLAG, Boolean.class, true);


    //Edit distance factory
    registerFlag(EDIT_DIST_LOGGING_AMOUNT_FLAG, Integer.class, 0);
    registerFlag(EDIT_DIST_HEURISTIC_ALIGNERS_FLAG, Boolean.class, true);
    registerFlag(EDIT_DIST_GOTOH_ONLY_FLAG);
    registerFlag(EDIT_DIST_SINGLE_INDEL_SEEDED_ONLY_FLAG);
    registerFlag(EDIT_DIST_INDEL_TABLE_FLAG, String.class, "");
    registerFlag(EDIT_DIST_LOG_AS_HISTOGRAM_FLAG);


    registerFlag(TEMP_FILES_DUMP_ALIGN_STATS_FLAG);
    registerFlag(MAP_KEEP_TEMPORARY_FILES);
    registerFlag(SLIDING_WINDOW_MAX_HITS_PER_POS_FLAG, Integer.class, 1000);

    //SAM
    registerFlag(SAM_ALLOW_FALLBACK_FOR_NON_INDEXED_REGIONS);
    registerFlag(GZIP_LEVEL, Integer.class, AdjustableGZIPOutputStream.DEFAULT_GZIP_LEVEL);

    // Aligners / all-paths
    registerFlag(DEFAULT_INDEL_LENGTH_FLAG, Integer.class, 7);

    registerFlag(ASSEMBLER_MAX_HITS_PER_START_POS_FLAG, Integer.class, 5);
    registerFlag(ASSEMBLER_INSERT_DEVIATIONS_FLAG, Integer.class, 4);


    // Complex caller
    registerFlag(COMPLEX_CALLER_MAX_HYPOTH_FLAG, Integer.class, 20);
    registerFlag(COMPLEX_REGION_INDEL_EXTENSION);
    registerFlag(COMPLEX_REGION_SIMPLE_REPEAT_LIMIT, Integer.class, 3);
    registerFlag(COMPLEX_REGION_SIMPLE_REPEAT_IMPL, String.class, "default");
    registerFlag(DUMP_COMPLEX_TRIGGER_SIGNALS, Boolean.class, false);
    registerFlag(SOFT_CLIP_COMPLEX_TRIGGER, Boolean.class, false);
    registerFlag(COMPLEX_EVIDENCE_DETAILS);
    registerFlag(COMPLEX_HYPOTHESIS_DETAILS);
    registerFlag(COMPLEX_HYPOTHESIS_NEW_PRIORS, Boolean.class, true);
    registerFlag(COMPLEX_HYPOTHESIS_ADJUST_PRIORS, Boolean.class, true);
    registerFlag(COMPLEX_CALLER_UNROLL_CG_FLAG, Boolean.class, true);

    // Misc calling
    registerFlag(CALLER_N_MIN_DEPTH, Integer.class, 5);
    registerFlag(FAMILY_CALLER_FALLBACK_FLAG, Boolean.class, false);
    registerFlag(FAMILY_CALLER_SORTED_HYPOTHESES, Boolean.class, true);
    registerFlag(MIN_BASE_QUALITY_AS_TWO, Boolean.class, false);
    registerFlag(QUALITY_CALIBRATION_MIN_EVIDENCE, Integer.class, 10);
    registerFlag(QUALITY_CALIBRATION_COVARIATE_INTERSECTION, Boolean.class, false);
    registerFlag(TUMOR_ALLELE_BALANCE, String.class, "binomial");


    // AVR, training on missing instances increases time and experience indicates is a bad idea
    // when there are lots of missing values.
    registerFlag(AVR_TRAIN_ON_MISSING_VALUES, Boolean.class, false);
    registerFlag(AVR_ALLOW_UNDECLARED_ATTRIBUTES, Boolean.class, false);

    //Aview
    registerFlag(CP437_DELETES, Boolean.class, false);
  }
}
