/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

/**
 * CLI class for <code>SamValidator</code>.
 */
public final class SamValidatorCli extends AbstractCli {

  private static final String MODULE_NAME = "samstats";
  private static final String PER_FILE_STATS = "per-file";
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
    flags.setDescription("Prints alignment statistics from the contents of the output SAM/BAM file.");
    CommonFlagCategories.setCategories(flags);
    CommonFlags.initReferenceTemplate(flags, true);
    final Flag<File> input = flags.registerRequired(File.class, CommonFlags.FILE, "SAM/BAM result file (must contain read-ids not read names)").setCategory(INPUT_OUTPUT);
    input.setMinCount(0).setMaxCount(Integer.MAX_VALUE);
    final Flag<File> inputList = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, CommonFlags.FILE, "file containing a list of SAM/BAM format files (1 per line) containing mapped reads").setCategory(INPUT_OUTPUT);
    flags.registerOptional('r', READS_DIR, File.class, CommonFlags.SDF, "reads SDF").setCategory(INPUT_OUTPUT);
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
      if (!CommonFlags.validateSDF(flags, CommonFlags.TEMPLATE_FLAG)
        || (flags.isSet(READS_DIR) && !CommonFlags.validateSDF(flags, READS_DIR))) {
        return false;
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
    final File template = (File) flags.getValue(CommonFlags.TEMPLATE_FLAG);
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
