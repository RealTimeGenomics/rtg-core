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

package com.rtg.assembler;

import static com.rtg.assembler.DeBruijnAssemblerCli.DIPLOID_RATIO;
import static com.rtg.assembler.DeBruijnAssemblerCli.KMER_SIZE;
import static com.rtg.assembler.GraphMapCli.MAX_INSERT;
import static com.rtg.assembler.GraphMapCli.MIN_INSERT;
import static com.rtg.assembler.GraphMapCli.MISMATCHES;
import static com.rtg.launcher.CommonFlags.END_READ_ID;
import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.launcher.CommonFlags.START_READ_ID;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.List;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ParamsCli;
import com.rtg.ngs.MapFlags;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.IORunnable;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;


/**
 */
public class AssembleCli extends ParamsCli<AssembleParams> {

  @Override
  public String moduleName() {
    return "assemble";
  }

  @Override
  public String description() {
    return "assemble reads into long sequences";
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(OUTPUT_FLAG);
  }

  @Override
  protected IORunnable task(AssembleParams params, OutputStream out) {
    return new AssembleTask(params, out, new GraphMapStatistics(params.directory()));
  }

  @Override
  protected AssembleParams makeParams() throws InvalidParamsException, IOException {
    return makeParamsLocal(mFlags);
  }

  @Override
  protected void initFlags() {
    initLocalFlags(mFlags);
  }

  protected static void initLocalFlags(CFlags flags) {
    CommonFlagCategories.setCategories(flags);
    flags.setValidator(new AssembleFlagsValidator());
    DeBruijnAssemblerCli.initCommonFlags(flags);
    GraphMapCli.initCommonFlags(flags);
    FilterPaths.initCommonFlags(flags);
    Consensus.initCommonFlags(flags);
    flags.registerOptional(GraphMapCli.ALIGNMENTS, "output an alignment file").setCategory(INPUT_OUTPUT);
    flags.registerOptional('g', GraphMapCli.GRAPH_FLAG, File.class, "Dir", "graph of the assembly to map against").setCategory(INPUT_OUTPUT);
  }

  private static class AssembleFlagsValidator implements Validator {

    private long sequences(final File dir) throws IOException {
      try (SequencesReader r = SequencesReaderFactory.createDefaultSequencesReader(dir)) {
        return r.numberSequences();
      }
    }

    @Override
    public boolean isValid(CFlags flags) {
      if (!CommonFlags.validateOutputDirectory(flags)) {
        return false;
      }
      try {
        final List<File> files = CommonFlags.getFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, true);
        files.addAll(CommonFlags.getFileList(flags, GraphMapCli.INPUT_LIST_FLAG_MATE_PAIR, GraphMapCli.MATE_PAIR, true));
        files.addAll(CommonFlags.getFileList(flags, GraphMapCli.INPUT_LIST_FLAG_454, GraphMapCli.FOUR_FIVE_FOUR, true));
        if (flags.isSet(START_READ_ID) || flags.isSet(END_READ_ID)) {
          if (files.size() != 1) {
            flags.setParseMessage("Can only specify read range with a single input set of reads");
            return false;
          }
        }
        if (!flags.isSet(GraphMapCli.GRAPH_FLAG)) {
          if (files.size() == 0) {
            flags.setParseMessage("You must supply at least one input read set or a pre-built graph");
            return false;
          }
          long totalSequences = 0;
          for (final File f : files) {
            if (ReaderUtils.isPairedEndDirectory(f)) {
              totalSequences += sequences(ReaderUtils.getLeftEnd(f));
              totalSequences += sequences(ReaderUtils.getRightEnd(f));
            } else {
              totalSequences += sequences(f);
            }
          }
          if (totalSequences == 0) {
            flags.setParseMessage("No reads in supplied input files");
            return false;
          }
        }
      } catch (IOException e) {
        flags.setParseMessage("Couldn't read the file list");
        return false;
      }
      if (!flags.checkInRange(MapFlags.WORDSIZE_FLAG, 1, Integer.MAX_VALUE)) {
        return false;
      }
      if (!flags.checkInRange(MapFlags.STEP_FLAG, 1, Integer.MAX_VALUE)) {
        return false;
      }
      if (!flags.checkInRange(KMER_SIZE, 1, Integer.MAX_VALUE)) {
        return false;
      }
      final IntegerOrPercentage maxMatedScore = (IntegerOrPercentage) flags.getValue(MISMATCHES);
      if (maxMatedScore.getValue(100) < 0) {
        Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, "--" + MISMATCHES, maxMatedScore + "", "0");
        return false;
      }
      if (flags.isSet(MAX_INSERT) ^ flags.isSet(MIN_INSERT)) {
        flags.setParseMessage("You must supply both --" + MAX_INSERT + " and --" + MIN_INSERT + " if you supply one");
        return false;
      }
      if (flags.isSet(MAX_INSERT) && (Integer) flags.getValue(MAX_INSERT) < (Integer) flags.getValue(MIN_INSERT)) {
        flags.setParseMessage("--" + MAX_INSERT + " should be larger than --" + MIN_INSERT);
        return false;
      }
      if ((Integer) flags.getValue(KMER_SIZE) < 1) {
        flags.setParseMessage("--" + KMER_SIZE + " should be positive");
        return false;
      }
      if (flags.isSet(DIPLOID_RATIO)) {
        final double ratio = (Double) flags.getValue(DIPLOID_RATIO);
        if (ratio > 1 || ratio < 0) {
          flags.setParseMessage("--" + DIPLOID_RATIO + " should be between 0 and 1");
        }
      }
      if (!CommonFlags.validateStartEnd(flags, START_READ_ID, END_READ_ID)) {
        return false;
      }
      return CommonFlags.validateThreads(flags);
    }
  }


  protected static AssembleParams makeParamsLocal(CFlags flags) throws IOException {
    final AssembleParams.Builder builder = AssembleParams.builder();
    builder.directory((File) flags.getValue(OUTPUT_FLAG))
        .consensusThreshold((Integer) flags.getValue(Consensus.CONSENSUS_READS))
        .graph((File) flags.getValue(GraphMapCli.GRAPH_FLAG))
        .reads(CommonFlags.getFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, true))
        .reads454(CommonFlags.getFileList(flags, GraphMapCli.INPUT_LIST_FLAG_454, GraphMapCli.FOUR_FIVE_FOUR, true))
        .readsMatePair(CommonFlags.getFileList(flags, GraphMapCli.INPUT_LIST_FLAG_MATE_PAIR, GraphMapCli.MATE_PAIR, true))
        .wordSize((Integer) flags.getValue(MapFlags.WORDSIZE_FLAG))
        .stepSize((Integer) flags.getValue(MapFlags.STEP_FLAG))
        .maxMismatches((IntegerOrPercentage) flags.getValue(MISMATCHES))
        .kmerSize((Integer) flags.getValue(KMER_SIZE))
        .minHashFrequency((Integer) flags.getValue(DeBruijnAssemblerCli.MIN_HASH_FREQUENCY))
        .mergeRatio((Double) flags.getValue(DIPLOID_RATIO))
        .minPathReads(flags.isSet(FilterPaths.MIN_PATH) ? (Integer) flags.getValue(FilterPaths.MIN_PATH) : -1)
        .minReadCount(flags.isSet(FilterPaths.READ_COUNT) ? (Integer) flags.getValue(FilterPaths.READ_COUNT) : -1)
        .alignments(flags.isSet(GraphMapCli.ALIGNMENTS))
        .minInsertSize(flags.isSet(MIN_INSERT) ? (Integer) flags.getValue(MIN_INSERT) : Integer.MAX_VALUE)
        .maxInsertSize(flags.isSet(MAX_INSERT) ? (Integer) flags.getValue(MAX_INSERT) : Integer.MIN_VALUE)
        .region(CommonFlags.getReaderRestriction(flags))
        .numberThreads(CommonFlags.parseThreads((Integer) flags.getValue(CommonFlags.THREADS_FLAG)));

    return builder.create();

  }

}
