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

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.GraphFactory;
import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.io.GraphReader;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ParamsCli;
import com.rtg.ngs.MapFlags;
import com.rtg.util.IORunnable;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.store.StoreDirProxy;
import com.rtg.util.store.StoreDirectory;

/**
 */
public class GraphMapCli extends ParamsCli<GraphMapParams> {

  static final String GRAPH_FLAG = "graph";
  static final String MISMATCHES = "mismatches";
  static final String MIN_INSERT = "min-insert";
  static final String MAX_INSERT = "max-insert";
  static final String MODULE_NAME = "graphmap";
  static final String FOUR_FIVE_FOUR = "454";
  static final String MATE_PAIR = "mate-pair";
  static final String INPUT_LIST_FLAG_454 = "input-list-454";
  static final String INPUT_LIST_FLAG_MATE_PAIR = "input-list-mate-pair";
  static final String ALIGNMENTS = "Xalignments";

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return null;
  }

  @Override
  protected IORunnable task(GraphMapParams params, OutputStream out) {
    return new GraphMapTask(params, out);
  }

  @Override
  protected GraphMapParams makeParams() throws IOException {
    return makeParamsLocal(mFlags);
  }

  protected static GraphMapParams makeParamsLocal(CFlags flags) throws IOException {
    final GraphMapParams.Builder builder = GraphMapParams.builder();
    builder.directory((File) flags.getValue(CommonFlags.OUTPUT_FLAG))
        .reads(CommonFlags.getFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, true))
        .reads454(CommonFlags.getFileList(flags, INPUT_LIST_FLAG_454, FOUR_FIVE_FOUR, true))
        .readsMatePair(CommonFlags.getFileList(flags, INPUT_LIST_FLAG_MATE_PAIR, MATE_PAIR, true))
        .wordSize((Integer) flags.getValue(MapFlags.WORDSIZE_FLAG))
        .stepSize((Integer) flags.getValue(MapFlags.STEP_FLAG))
        .graph(loadGraph(new StoreDirProxy((File) flags.getValue(GRAPH_FLAG))))
        .numberThreads(CommonFlags.parseThreads((Integer) flags.getValue(CommonFlags.THREADS_FLAG)))
        .maxMismatches((IntegerOrPercentage) flags.getValue(MISMATCHES));
    if (flags.isSet(MIN_INSERT) && flags.isSet(MAX_INSERT)) {
      builder
        .minInsertSize((Integer) flags.getValue(MIN_INSERT))
        .maxInsertSize((Integer) flags.getValue(MAX_INSERT));
    }
    if (flags.isSet(ALIGNMENTS)) {
      builder.alignmentFile((File) flags.getValue(ALIGNMENTS));
    }
     return builder.create();
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  @Override
  protected void initFlags() {
    initLocalFlags(mFlags);
  }

  protected static void initCommonFlags(CFlags flags) {
    CommonFlags.initThreadsFlag(flags);
    flags.registerOptional('w', MapFlags.WORDSIZE_FLAG, Integer.class, CommonFlags.INT, "word size", 18).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('s', MapFlags.STEP_FLAG, Integer.class, CommonFlags.INT, "step size", 18).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('a', MISMATCHES, IntegerOrPercentage.class, CommonFlags.INT, "number of bases that may mismatch in an alignment or percentage of read that may mismatch", new IntegerOrPercentage(0)).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('m', MIN_INSERT, Integer.class, CommonFlags.INT, "minimum insert size between fragments").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('M', MAX_INSERT, Integer.class, CommonFlags.INT, "maximum insert size between fragments").setCategory(SENSITIVITY_TUNING);
    final Flag<File> reads454 = flags.registerOptional('f', FOUR_FIVE_FOUR, File.class, CommonFlags.SDF, "SDF containing 454 reads").setCategory(INPUT_OUTPUT);
    final Flag<File> listFlag454 = flags.registerOptional('F', INPUT_LIST_FLAG_454, File.class, CommonFlags.FILE, "file containing a list of SDF directories (1 per line) containing 454 sequences to assemble").setCategory(INPUT_OUTPUT);
    reads454.setMinCount(0);
    reads454.setMaxCount(Integer.MAX_VALUE);
    final Flag<File> readsMatePair = flags.registerOptional('j', MATE_PAIR, File.class, CommonFlags.SDF, "SDF containing mate pair reads").setCategory(INPUT_OUTPUT);
    final Flag<File> listFlagMatePair = flags.registerOptional('J', INPUT_LIST_FLAG_MATE_PAIR, File.class, CommonFlags.FILE, "file containing a list of SDF directories (1 per line) containing mate pair sequences to assemble").setCategory(INPUT_OUTPUT);
    readsMatePair.setMinCount(0);
    readsMatePair.setMaxCount(Integer.MAX_VALUE);
    CommonFlagCategories.setCategories(flags);
    final Flag<File> inFlag = flags.registerRequired(File.class, CommonFlags.SDF, "SDF directories containing reads to map");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    final Flag<File> listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, CommonFlags.FILE, "file containing a list of SDF directories (1 per line) containing sequences to assemble").setCategory(INPUT_OUTPUT);
    flags.addRequiredSet(inFlag);
    flags.addRequiredSet(listFlag);
    flags.addRequiredSet(reads454);
    flags.addRequiredSet(listFlag454);
    flags.addRequiredSet(readsMatePair);
    flags.addRequiredSet(listFlagMatePair);
  }

  protected static void initLocalFlags(CFlags flags) {
    flags.setValidator(new GraphMapValidator());
    initCommonFlags(flags);
    CommonFlags.initOutputDirFlag(flags);
    flags.registerRequired('g', GRAPH_FLAG, File.class, "Dir", "graph of the assembly to map against").setCategory(INPUT_OUTPUT);
    flags.registerOptional(ALIGNMENTS, File.class, CommonFlags.FILE, "alignments will be written to this file").setCategory(INPUT_OUTPUT);
  }

  /**
   * Noddy main.
   * @param args see usage
   */
  public static void main(final String[] args) {
    new GraphMapCli().mainInit(args, System.out, System.err);
  }

  private static class GraphMapValidator implements Validator {
    /**
     * Check the file list and anonymous file input flags.
     * @param flags the flags to check
     * @return <code>true</code> if all okay <code>false</code> otherwise
     */
    public static boolean checkSdfFileList(CFlags flags) {
      final Collection<File> files;
      try {
        files = CommonFlags.getFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, true);
      } catch (final IOException e) {
        flags.setParseMessage("An error occurred reading " + flags.getValue(CommonFlags.INPUT_LIST_FLAG));
        return false;
      }
      try {
        files.addAll(CommonFlags.getFileList(flags, INPUT_LIST_FLAG_454, FOUR_FIVE_FOUR, true));
      } catch (final IOException e) {
        flags.setParseMessage("An error occurred reading " + flags.getValue(INPUT_LIST_FLAG_454));
        return false;
      }
      if (files.size() == 0) {
        flags.setParseMessage("No input files specified.");
        return false;
      }
      return true;
    }

    @Override
    public boolean isValid(CFlags flags) {
      if (!CommonFlags.validateOutputDirectory(flags)) {
        return false;
      }
      if (!checkSdfFileList(flags)) {
        return false;
      }
      if (flags.isSet(ALIGNMENTS) && ((File) flags.getValue(ALIGNMENTS)).exists()) {
        flags.setParseMessage("Specified alignment file already exists");
        return false;
      }
      if (!flags.checkInRange(MapFlags.WORDSIZE_FLAG, 1, Integer.MAX_VALUE)) {
        return false;
      }
      if (!flags.checkInRange(MapFlags.STEP_FLAG, 1, Integer.MAX_VALUE)) {
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
      if (flags.isSet(MAX_INSERT) && flags.isSet(MIN_INSERT) && (Integer) flags.getValue(MAX_INSERT) < (Integer) flags.getValue(MIN_INSERT)) {
          flags.setParseMessage("--" + MAX_INSERT + " should be larger than --" + MIN_INSERT);
          return false;
      }

      return CommonFlags.validateThreads(flags);
    }
  }

  private static MutableGraph loadGraph(StoreDirectory graphFile) throws IOException {
    final Map<String, String> contigAttrs = new HashMap<>();
    final Map<String, String> pathAttrs = new HashMap<>();
    contigAttrs.put(GraphKmerAttribute.READ_COUNT, "number of reads that map uniquely within this contig");
    pathAttrs.put(GraphKmerAttribute.READ_COUNT, "number of reads that map uniquely across the links in this path");
    final Graph graph = GraphReader.read(GraphFactory.KMER, graphFile, contigAttrs, pathAttrs);
    assert graph instanceof MutableGraph;
    return (MutableGraph) graph;
  }

}
