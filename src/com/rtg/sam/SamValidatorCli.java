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

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Collection;

import com.rtg.alignment.EditDistanceFactory;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.reader.ReaderUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;

/**
 * CLI class for <code>SamValidator</code>.
 */
public final class SamValidatorCli extends AbstractCli {

  private static final String MODULE_NAME = "samstats";
  private static final String PER_FILE_STATS = "per-file";
  private static final String TEMPLATE_DIR = "template";
  private static final String READS_DIR = "reads";
  private static final String VALIDATE_FLAG = "validate";
  private static final String CONSENSUS_FLAG = "consensus";
  private static final String DISTRIBUTIONS_FLAG = "distributions";
  private static final String IGNORE_CG_FRAGMENT_LENGTH = "Xignore-cg-fragment";

  static final String GAP_OPEN_PENALTY_FLAG = "Xgap-open-penalty";
  static final String GAP_EXTEND_PENALTY_FLAG = "Xgap-extend-penalty";
  static final String MISMATCH_PENALTY_FLAG = "Xmismatch-penalty";
  static final String UNKNOWNS_PENALTY_FLAG = "Xunknowns-penalty";

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return "print statistics about a SAM/BAM file";
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }
  /**
   * construct a flags object
   * @param flags flags object to set up
   */
  protected void initFlags(final CFlags flags) {
    flags.registerExtendedHelp();
    flags.setDescription("Prints alignment statistics from the contents of the output SAM/BAM file.");
    CommonFlagCategories.setCategories(flags);
    flags.registerRequired('t', TEMPLATE_DIR, File.class, "SDF", "template SDF").setCategory(INPUT_OUTPUT);
    final Flag input = flags.registerRequired(File.class, "FILE", "SAM/BAM result file (must contain read-ids not read names)").setCategory(INPUT_OUTPUT);
    input.setMinCount(0).setMaxCount(Integer.MAX_VALUE);
    final Flag inputList = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, "FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);
    flags.registerOptional('r', READS_DIR, File.class, "SDF", "reads SDF").setCategory(INPUT_OUTPUT);
    flags.registerOptional(VALIDATE_FLAG, "validate mapping of read to reference. Tests matching of bases according to CIGAR format").setCategory(REPORTING);
    flags.registerOptional(CONSENSUS_FLAG, "record consensus data. Requires roughly 5 times template length of RAM").setCategory(REPORTING);
    flags.registerOptional('D', DISTRIBUTIONS_FLAG, "display distributions of insert sizes, alignment scores and read hits").setCategory(REPORTING);
    flags.registerOptional(IGNORE_CG_FRAGMENT_LENGTH, "ignore unusual Complete Genomics fragment lengths.").setCategory(REPORTING);
    flags.registerOptional(PER_FILE_STATS, "output per-file statistics").setCategory(REPORTING);

    flags.registerOptional(GAP_OPEN_PENALTY_FLAG, Integer.class, CommonFlags.INT, "penalty for a gap open during alignment score checking").setCategory(REPORTING);
    flags.registerOptional(GAP_EXTEND_PENALTY_FLAG, Integer.class, CommonFlags.INT, "penalty for a gap extension during alignment score checking").setCategory(REPORTING);
    flags.registerOptional(MISMATCH_PENALTY_FLAG, Integer.class, CommonFlags.INT, "penalty for a mismatch during alignment score checking").setCategory(REPORTING);
    flags.registerOptional(UNKNOWNS_PENALTY_FLAG, Integer.class, CommonFlags.INT, "penalty for an unknown nucleotide during alignment score checking").setCategory(REPORTING);

    flags.addRequiredSet(input);
    flags.addRequiredSet(inputList);
    flags.setValidator(VALIDATOR);
  }

  private static final Validator VALIDATOR = new Validator() {
    @Override
    public boolean isValid(final CFlags flags) {
      if (!CommonFlags.checkFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, Integer.MAX_VALUE)) {
        return false;
      }
      final String templateFile = flags.getFlag(TEMPLATE_DIR).getValue().toString();
      final File template = new File(templateFile);
      if (!template.exists()) {
        Diagnostic.error(ErrorType.INFO_ERROR, "The specified SDF, \"" + template.getPath() + "\", does not exist.");
        return false;
      }
      if (!template.isDirectory()) {
        Diagnostic.error(ErrorType.INFO_ERROR, "The specified file, \"" + template.getPath() + "\", is not an SDF.");
        return false;
      }
      if (flags.getFlag(READS_DIR).isSet()) {
        final String readsFile = flags.getFlag(READS_DIR).getValue().toString();
        final File reads = new File(readsFile);
        if (!reads.exists()) {
          Diagnostic.error(ErrorType.INFO_ERROR, "The specified SDF, \"" + reads.getPath() + "\", does not exist.");
          return false;
        }
        if (!reads.isDirectory()) {
          Diagnostic.error(ErrorType.INFO_ERROR, "The specified file, \"" + reads.getPath() + "\", is not an SDF.");
          return false;
        }
      }
      if (flags.isSet(GAP_EXTEND_PENALTY_FLAG) || flags.isSet(GAP_OPEN_PENALTY_FLAG) || flags.isSet(MISMATCH_PENALTY_FLAG) || flags.isSet(UNKNOWNS_PENALTY_FLAG)) {
        if (!(flags.isSet(GAP_EXTEND_PENALTY_FLAG) && flags.isSet(GAP_OPEN_PENALTY_FLAG) && flags.isSet(MISMATCH_PENALTY_FLAG) && flags.isSet(UNKNOWNS_PENALTY_FLAG))) {
          Diagnostic.error("Must specify all of --" + MISMATCH_PENALTY_FLAG + ", --" + GAP_OPEN_PENALTY_FLAG + ", --" + GAP_EXTEND_PENALTY_FLAG + ", --" + UNKNOWNS_PENALTY_FLAG + " if specifying any.");
          return false;
        }
      }
      return true;

    }
  };

  @Override
  protected int mainExec(final OutputStream out, final PrintStream err) throws IOException {
    final CFlags flags = mFlags;

    final Collection<File> inputFiles = CommonFlags.getFileList(mFlags, CommonFlags.INPUT_LIST_FLAG, null, false);
    final File template = (File) flags.getValue(TEMPLATE_DIR);
    final boolean validate = flags.isSet(VALIDATE_FLAG);
    final boolean showConsensus = flags.isSet(CONSENSUS_FLAG);
    final boolean printHistograms = flags.isSet(DISTRIBUTIONS_FLAG);
    final boolean ignoreCgFragmentSize = flags.isSet(IGNORE_CG_FRAGMENT_LENGTH);
    final boolean perFileStats = flags.isSet(PER_FILE_STATS);

    final PrintStream outps = new PrintStream(out);
    try {
      final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(flags.isSet(GAP_OPEN_PENALTY_FLAG) ? (Integer) flags.getValue(GAP_OPEN_PENALTY_FLAG) : EditDistanceFactory.DEFAULT_GAP_OPEN_PENALTY)
                                                     .gapExtendPenalty(flags.isSet(GAP_EXTEND_PENALTY_FLAG) ? (Integer) flags.getValue(GAP_EXTEND_PENALTY_FLAG) : EditDistanceFactory.DEFAULT_GAP_EXTEND_PENALTY)
                                                     .substitutionPenalty(flags.isSet(MISMATCH_PENALTY_FLAG) ? (Integer) flags.getValue(MISMATCH_PENALTY_FLAG) : EditDistanceFactory.DEFAULT_SUBSTITUTION_PENALTY)
                                                     .unknownsPenalty(flags.isSet(UNKNOWNS_PENALTY_FLAG) ? (Integer) flags.getValue(UNKNOWNS_PENALTY_FLAG) : EditDistanceFactory.DEFAULT_UNKNOWNS_PENALTY)
                                                     .create();
      final SamValidator sv = new SamValidator(outps, err, validate, showConsensus, printHistograms, ignoreCgFragmentSize, perFileStats, params, flags.isSet(MISMATCH_PENALTY_FLAG)); //if any penalties flags are set, all must be, so only checking mismatch is ok

      final File reads = (File) flags.getValue(READS_DIR);
      final File left;
      final File right;
      if (ReaderUtils.isPairedEndDirectory(reads)) {
        left = ReaderUtils.getLeftEnd(reads);
        right = ReaderUtils.getRightEnd(reads);
      } else {
        left = reads;
        right = null;
      }
      sv.checkSAMAlign(template, inputFiles, left, right);
    } finally {
      outps.flush();
    }

    Diagnostic.closeLog();
    Diagnostic.deleteLog();
    return 0;
  }


  /**
   * Main program for getting sam stats
   * @param args command line arguments.
   */
  public static void main(final String[] args) {
    new SamValidatorCli().mainExit(args);
  }
}
