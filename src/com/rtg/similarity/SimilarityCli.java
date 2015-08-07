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
package com.rtg.similarity;

import static com.rtg.launcher.BuildCommon.RESOURCE;
import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.ngs.MapFlags.STEP_FLAG;
import static com.rtg.ngs.MapFlags.WORDSIZE_FLAG;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.index.Index;
import com.rtg.index.IndexUtils;
import com.rtg.index.hash.ExactHashFunction;
import com.rtg.index.hash.HashLoop;
import com.rtg.index.hash.IncrementalHashLoop;
import com.rtg.index.params.CountParams;
import com.rtg.index.params.ParamsUtils;
import com.rtg.index.similarity.BinaryTree;
import com.rtg.index.similarity.IndexSimilarity;
import com.rtg.index.similarity.NeighborJoining;
import com.rtg.index.similarity.SimilarityMatrix;
import com.rtg.launcher.BuildParams;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.GlobalFlags;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.NoStatistics;
import com.rtg.launcher.ParamsCli;
import com.rtg.launcher.ParamsTask;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.ProgramMode;
import com.rtg.ngs.MapFlags;
import com.rtg.reader.IndexFile;
import com.rtg.reader.PrereadNamesInterface;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.similarity.BuildSearchParams.BuildSearchParamsBuilder;
import com.rtg.taxonomy.TaxonNode;
import com.rtg.taxonomy.Taxonomy;
import com.rtg.taxonomy.TaxonomyUtils;
import com.rtg.usage.UsageMetric;
import com.rtg.util.IORunnable;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.MathUtils;
import com.rtg.util.Pair;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.OneShotTimer;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.intervals.LongRange;

/**
 * The similarity command estimates how closely related a set of samples are to each other.
 * It produces output in the form of a similarity matrix, a principal component analysis
 * and two formats for phylogenetic trees showing relationships.
 */
public final class SimilarityCli extends ParamsCli<BuildSearchParams> {

  private static final String UNIQUE_WORDS = "unique-words";
  /** The file name of where the modified New Hampshire tree is written. */
  static final String TREE_SUFFIX = "closest.tre";
  static final String XML_SUFFIX = "closest.xml";
  /** The file name of where the similarity matrix is written. */
  public static final String SIMILARITY_SUFFIX = "similarity.tsv";
  /** The file name of where the similarity PCA output is written. */
  public static final String PCA_SUFFIX = "similarity.pca";

  static final String MODULE_NAME = "similarity";
  private static final String INPUT_FLAG = "input";

  private static final String MAX_READS_FLAG = "max-reads";

  static final int DEFAULT_WORD_SIZE = 20;
  static final int DEFAULT_STEP_SIZE = 1;

  private static final Validator VALIDATOR = new SimilarityFlagsValidator();

  @TestClass(value = {"com.rtg.similarity.SimilarityFlagsValidatorTest"})
  private static class SimilarityFlagsValidator implements Validator {
    @Override
    public boolean isValid(final CFlags flags) {
      if (flags.isSet(INPUT_FLAG)) {
        if (flags.isSet(CommonFlags.INPUT_LIST_FLAG)) {
          flags.setParseMessage("Only set one of --" + INPUT_FLAG + " or --" + CommonFlags.INPUT_LIST_FLAG);
          return false;
        }
        if (flags.isSet(MAX_READS_FLAG)) {
          flags.setParseMessage("Only set --" + MAX_READS_FLAG + " when using --" + CommonFlags.INPUT_LIST_FLAG);
          return false;
        }
        if (!CommonFlags.validateSDF(flags, INPUT_FLAG)) {
          return false;
        }
        final File input = (File) flags.getValue(INPUT_FLAG);
        if (ReaderUtils.isPairedEndDirectory(input)) {
          Diagnostic.error(ErrorType.INFO_ERROR, "The specified SDF, \"" + input.getPath() + "\", is a paired end SDF.");
          return false;
        }
      } else if (!flags.isSet(CommonFlags.INPUT_LIST_FLAG)) {
        flags.setParseMessage("Must set one of --" + INPUT_FLAG + " or --" + CommonFlags.INPUT_LIST_FLAG);
        return false;
      } else {
        final File inputList = (File) flags.getValue(CommonFlags.INPUT_LIST_FLAG);
        if (!inputList.exists()) {
          flags.setParseMessage("The specified list file, \"" + inputList.getPath() + "\", does not exist.");
          return false;
        }
        if (inputList.isDirectory()) {
          flags.setParseMessage("The specified list file, \"" + inputList.getPath() + "\", is a directory.");
          return false;
        }
        if (flags.isSet(MAX_READS_FLAG) && (Integer) flags.getValue(MAX_READS_FLAG) < 1) {
          flags.setParseMessage("The --" + MAX_READS_FLAG + " must be greater than 0");
          return false;
        }
      }
      if (!CommonFlags.validateOutputDirectory(flags)) {
        return false;
      }
      if (flags.isSet(WORDSIZE_FLAG)) {
        if (!CommonFlags.validateFlagBetweenValues(flags, WORDSIZE_FLAG, 1, 32)) {
          return false;
        }
      }
      if (!MapFlags.validateStepAndWordSize(flags)) {
        return false;
      }
      return true;
    }
  }

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return "calculate similarity matrix and nearest neighbor tree";
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(OUTPUT_FLAG);
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  protected static void initFlags(CFlags flags) {
    flags.setValidator(VALIDATOR);
    flags.setDescription("Produces a similarity matrix and nearest neighbor tree from the input sequences or reads.");
    CommonFlagCategories.setCategories(flags);
    CommonFlags.initOutputDirFlag(flags);
    final Flag inFlag = flags.registerOptional('i', INPUT_FLAG, File.class, "SDF", RESOURCE.getString("SUBJECT_DESC")).setCategory(INPUT_OUTPUT);
    final Flag listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, "FILE", "file containing a labeled list of SDF files (1 label and file per line format:[label][space][file])").setCategory(INPUT_OUTPUT);
    MapFlags.initWordSize(flags, "word size (Default is " + DEFAULT_WORD_SIZE + ")");
    MapFlags.initStepSize(flags, "step size (Default is " + DEFAULT_STEP_SIZE + ")");
    flags.registerOptional(UNIQUE_WORDS, "count only unique words").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional(MAX_READS_FLAG, Integer.class, "INT", "maximum number of reads to use from each input SDF").setCategory(CommonFlagCategories.UTILITY);
    flags.addRequiredSet(inFlag);
    flags.addRequiredSet(listFlag);
  }

  protected static class SimilarityTask extends ParamsTask<BuildSearchParams, NoStatistics> {

    public SimilarityTask(final BuildSearchParams params, final OutputStream out, final UsageMetric usageMetric) {
      super(params, out, new NoStatistics(), usageMetric);
    }

    @Override
    protected void exec() {
      buildAndSearch(mParams);
    }
  }

  /**
   * Get  a string with a  human readable description of the memory usage.
   * It conforms to the usage that each line starts with the word "Memory" followed by
   * a unique identifying string that doesn't contain spaces followed by the
   * number of bytes. All entries are to be disjoint (so no totals are to be included).
   * @param buildSearchParams parameters to be used to create <code>BuildSearch</code>
   * @param bufferLength previously calculated buffer size.
   * @return a string with a  human readable description of the memory usage.
   */
  static String memToString(final BuildSearchParams buildSearchParams, final long bufferLength) {
    final StringBuilder sb = new StringBuilder();
    memToString(sb, buildSearchParams, bufferLength);
    return sb.toString();
  }

  /**
   * Get  a string with a  human readable description of the memory usage.
   * It conforms to the usage that each line starts with the word "Memory" followed by
   * a unique identifying string that doesn't contain spaces followed by the
   * number of bytes. All entries are to be disjoint (so no totals are to be included).
   * @param buildSearchParams parameters to be used to create <code>BuildSearch</code>
   * @param bufferLength previously calculated buffer size.
   * @param sb where to put the result.
   */
  static void memToString(final StringBuilder sb, final BuildSearchParams buildSearchParams, final long bufferLength) {
    sb.append(ParamsUtils.memToString("Shared_buffer", bufferLength));
    sb.append(IndexUtils.memString(buildSearchParams.build()));
  }

  private static class IncrementalIdMap extends HashMap<Integer, Integer> {
    private int mUnused = 0;
    int getId(final int v) {
      final Integer r = get(v);
      if (r != null) {
        return r;
      }
      put(v, mUnused);
      // System.err.println("Allocated " + mUnused + " to taxon id=" + v);
      return mUnused++;
    }
  }

  private static Pair<long[], List<String>> sdfToSimIdViaTaxId(final SequencesReader reader) throws IOException {
    if (TaxonomyUtils.hasTaxonomyInfo(reader)) {
      Diagnostic.userLog("Using taxonomy information for collating sequences");
      if (reader.numberSequences() > Integer.MAX_VALUE) {
        throw new UnsupportedOperationException("Too many sequences");
      }
      final long[] sdfIdToSimId = new long[(int) reader.numberSequences()];
      final List<String> taxonNames = new ArrayList<>();
      final IncrementalIdMap idMap = new IncrementalIdMap();
      final Map<String, Integer> seqNameToTaxonIdMap = TaxonomyUtils.loadTaxonomyMapping(reader);
      final Taxonomy taxonomy = TaxonomyUtils.loadTaxonomy(reader);
      final PrereadNamesInterface names = reader.names();
      long seen = -1;
      for (int k = 0; k < sdfIdToSimId.length; k++) {
        final int taxonId = seqNameToTaxonIdMap.get(names.name(k));
        sdfIdToSimId[k] = idMap.getId(taxonId);
        if (sdfIdToSimId[k] > seen) {
          // New taxon encountered
          seen = sdfIdToSimId[k];
          final TaxonNode tn = taxonomy.get(taxonId);
          taxonNames.add(tn.getName());
        }
      }
      Diagnostic.userLog("There are " + taxonNames.size() + " species in the reference set");
      return new Pair<>(sdfIdToSimId, taxonNames);
    } else {
      return null;
    }
  }

  /**
   * @param params parameters for the build and search.
   */
  public static void buildAndSearch(final BuildSearchParams params) {
    try {
      Diagnostic.userLog("Parameters:" + StringUtils.LS + params.toString());
      final long bufferLength = params.bufferLength();
      Diagnostic.userLog("Usage of memory" + StringUtils.LS + memToString(params, bufferLength));

      // Make all the components we need
      final IndexSimilarity index = new IndexSimilarity(params.build(), null, false, Integer.MAX_VALUE, 0, params.uniqueWords(), 1);

      // Search the queries and write hits
      final long numSequences;
      final byte[] buffer = makeBuffer(bufferLength);
      final List<String> names;
      if (params.build().sequences().directory() != null) {
        Diagnostic.progress("Input for index starting");
        try (final SequencesReader reader = params.build().sequences().reader()) {
          final Pair<long[], List<String>> sdfToSimId = sdfToSimIdViaTaxId(reader);
          final HashLoop buildLoop = makeBuild(index, params.build(), -1, sdfToSimId == null ? null : sdfToSimId.getA());
          buildLoop.execLoop(params.build().sequences(), buffer);
          index.freeze();
          Diagnostic.progress("Input for post-freeze index starting");
          buildLoop.execLoop(params.build().sequences(), buffer);
          if (sdfToSimId == null) {
            names = null;
            numSequences = params.build().sequences().numberSequences();
          } else {
            names = sdfToSimId.getB();
            numSequences = names.size();
          }
        }
      } else {
        for (int i = 0; i < params.sequences().size(); i++) {
          final Pair<String, List<SequenceParams>> pair = params.sequences().get(i);
          Diagnostic.progress("Input for index starting label \"" + pair.getA() +  "\" (" + (i + 1) + "/" + params.sequences().size() + ")");
          for (int j = 0; j < pair.getB().size(); j++) {
            Diagnostic.progress("Input for index label \"" + pair.getA() +  "\" starting directory (" + (j + 1) + "/" + pair.getB().size() + ")");
            final ISequenceParams seqParams = pair.getB().get(j);
            makeBuild(index, params.build(), i, null).execLoop(seqParams, buffer);
            seqParams.close();
          }
        }
        index.freeze();
        for (int i = 0; i < params.sequences().size(); i++) {
          final Pair<String, List<SequenceParams>> pair = params.sequences().get(i);
          Diagnostic.progress("Input for post-freeze index starting label \"" + pair.getA() +  "\" (" + (i + 1) + "/" + params.sequences().size() + ")");
          for (int j = 0; j < pair.getB().size(); j++) {
            Diagnostic.progress("Input for post-freeze index label \"" + pair.getA() +  "\" starting directory (" + (j + 1) + "/" + pair.getB().size() + ")");
            final ISequenceParams seqParams = pair.getB().get(j);
            makeBuild(index, params.build(), i, null).execLoop(seqParams, buffer);
            seqParams.close();
          }
        }
        names = null;
        numSequences = params.sequences().size();
      }

      // Build the internal indexes which allow searching to be done
      Diagnostic.progress("Input for index finished. Starting indexing");
      index.freeze();
      //System.err.println(index);
      Diagnostic.userLog("Memory performance " + StringUtils.LS + index.infoString());

      final boolean doPca = GlobalFlags.getBooleanValue(GlobalFlags.SIMILARITY_PCA_FLAG);

      try (final Writer writeSimi = new OutputStreamWriter(params.outStream(SIMILARITY_SUFFIX))) {
        try (final Writer writeTree = new OutputStreamWriter(params.outStream(TREE_SUFFIX))) {
          try (final Writer writeXml = new OutputStreamWriter(params.outStream(XML_SUFFIX))) {
            try (final Writer writePca = doPca ? new OutputStreamWriter(params.outStream(PCA_SUFFIX)) : null) {
              Diagnostic.progress("Indexing finished. Starting similarity.");
              similarity(index, numSequences, names, params, params.directory().toString(), writeSimi, writePca, writeTree, writeXml);
              Diagnostic.progress("Similarity finished.");
            }
          }
        }
      }
    } catch (final IOException e) {
      throw new SlimException(e);
    }
  }

  protected static HashLoop makeBuild(final Index index, final BuildParams buildParams, final int labelIndex, final long[] sdfIdToTaxonId) {
    final HashLoop subjectHashLoop;
    final boolean dualMode = true;
    final int winBits = buildParams.windowBits();
    if (winBits > 64) {
      throw new SlimException(ErrorType.INFO_ERROR, "Word size > 32");
    } else {
      final ExactHashFunction exf = new ExactHashFunction(buildParams, dualMode);
      if (sdfIdToTaxonId != null) {
        subjectHashLoop = new IncrementalHashLoop(buildParams.stepSize(), exf, dualMode) {
          @Override
          public void hashCall(final long hash, final int internalId, final int stepPosition) {
            throw new UnsupportedOperationException();
          }

          @Override
          public void hashCallBidirectional(final long hashForward, final long hashReverse, final int stepPosition, final int internalId) {
            if (hashForward < hashReverse) {
              index.add(hashForward, sdfIdToTaxonId[internalId]);
            } else {
              index.add(hashReverse, sdfIdToTaxonId[internalId]);
            }
          }

        };
      } else if (labelIndex < 0) {
        subjectHashLoop = new IncrementalHashLoop(buildParams.stepSize(), exf, dualMode) {
          @Override
          public void hashCall(final long hash, final int internalId, final int stepPosition) {
            throw new UnsupportedOperationException();
          }

          @Override
          public void hashCallBidirectional(final long hashForward, final long hashReverse, final int stepPosition, final int internalId) {
            //System.err.println("hashF=" + hashForward + " hashR=" + hashReverse + " id=" + internalId);
            if (hashForward < hashReverse) {
              index.add(hashForward, internalId);
            } else {
              index.add(hashReverse, internalId);
            }
          }
        };
      } else {
        subjectHashLoop = new IncrementalHashLoop(buildParams.stepSize(), exf, dualMode) {
          @Override
          public void hashCall(final long hash, final int internalId, final int stepPosition) {
            throw new UnsupportedOperationException();
          }

          @Override
          public void hashCallBidirectional(final long hashForward, final long hashReverse, final int stepPosition, final int internalId) {
            //System.err.println("hashF=" + hashForward + " hashR=" + hashReverse + " id=" + internalId + " labelIndex=" + labelIndex);
            if (hashForward < hashReverse) {
              index.add(hashForward, labelIndex);
            } else {
              index.add(hashReverse, labelIndex);
            }
          }

        };
      }
    }
    return subjectHashLoop;
  }

  /** Make a buffer long enough to be used for both build and search. */
  private static byte[] makeBuffer(final long maxSequence) {
    if (maxSequence > Integer.MAX_VALUE) {
      throw new SlimException(ErrorType.SEQUENCE_TOO_LONG, String.valueOf(maxSequence));
    }
    return new byte[(int) maxSequence];
  }

  static void similarity(final IndexSimilarity index, final long numSequences, List<String> presetNames, final BuildSearchParams params, final String outDir, final Appendable simiOut, final Appendable pcaOut, final Appendable treeOut, final Appendable xmlOut) throws IOException {
    final List<String> names;
    final OneShotTimer matrixTimer = new OneShotTimer("Ph_similarity_matrix");
    final SimilarityMatrix matrix = index.similarity(numSequences);
    //System.err.println(matrix);
    matrixTimer.stopLog();

    final OneShotTimer neighTimer = new OneShotTimer("Ph_similarity_neighbor");
    // 42 so reproduces the old versions behavior for regression
    final NeighborJoining neigh = new NeighborJoining(42);
    final BinaryTree tree;
    if (presetNames == null) {
      if (params.build().sequences().directory() != null) {
        names = makeNames(params.build().sequences().reader());
      } else {
        names = new ArrayList<>();
        for (int i = 0; i < params.sequences().size(); i++) {
          names.add(params.sequences().get(i).getA());
        }
      }
    } else {
      names = presetNames;
    }
    tree = neigh.neighborJoin(names, matrix);
    neighTimer.stopLog();

    try {
      if (treeOut != null) {
        tree.newick(treeOut);
      }
      if (xmlOut != null) {
        tree.phyloXml(xmlOut);
      }
    } catch (final IOException e) {
      throw new SlimException(e, ErrorType.WRITING_ERROR, outDir);
    }

    final OneShotTimer outTimer = new OneShotTimer("Ph_similarity_out");
    try {
      simiOut.append("#rtg ").append(CommandLine.getCommandLine());
      simiOut.append(StringUtils.LS);
      simiOut.append(matrix.toString(names));
    } catch (final IOException e) {
      throw new SlimException(e, ErrorType.WRITING_ERROR, outDir);
    }
    outTimer.stopLog();

    if (pcaOut != null) {
      // Do principle component analysis
      final OneShotTimer svdTimer = new OneShotTimer("Ph_SVD_out");
      final SimilaritySvd simMatrix = new SimilaritySvd(matrix.length());
      for (int j = 0; j < matrix.length(); j++) {
        for (int i = 0; i < matrix.length(); i++) {
          simMatrix.put(i, j, (long) matrix.get(i, j)); // matrix value are integers?
        }
        simMatrix.putName(j, names.get(j));
      }
      simMatrix.decompose(3);  // want 3 dimensions, mainly for plotting
      svdTimer.stopLog();

      try {
        for (int j = 0; j < simMatrix.getSvdRowsLength(); j++) {
          for (int i = 0; i < simMatrix.getSvdDimension(); i++) {
            pcaOut.append(Utils.realFormat(simMatrix.getSvdValue(j, i), 4)).append(StringUtils.TAB);
          }
          pcaOut.append(simMatrix.getSvdName(j)).append(StringUtils.LS);
        }
      } catch (final IOException e) {
        throw new SlimException(e, ErrorType.WRITING_ERROR, outDir);
      }
    }
  }

  /**
   * Get all sequence names in the reader.
   * @param reader sequences from where to get names.
   * @return list of all sequence names in the reader.
   * @throws IOException When IO errors occur
   */
  static ArrayList<String> makeNames(final SequencesReader reader) throws IOException {
    final ArrayList<String> nodeNames;
    nodeNames = new ArrayList<>();
    for (int i = 0; i < reader.numberSequences(); i++) {
      nodeNames.add(reader.name(i));
    }
    return nodeNames;
  }

  /**
   * Main program for building and searching. Use -h to get help.
   * @param args command line arguments.
   */
  public static void main(final String[] args) {
    new SimilarityCli().mainExit(args);
  }

  @Override
  protected BuildSearchParams makeParams() throws InvalidParamsException, IOException {
    final ProgramMode pm = ProgramMode.PHYLOGENY;
    final File output = (File) mFlags.getValue(OUTPUT_FLAG);
    final Integer window = (Integer) mFlags.getValue(WORDSIZE_FLAG);
    final Integer step = (Integer) mFlags.getValue(STEP_FLAG);

    final int nwindow = window == null ? DEFAULT_WORD_SIZE : window;
    final int nstep = step == null ? DEFAULT_STEP_SIZE : step;

    if (nstep > nwindow) {
      throw new InvalidParamsException(ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, "--" + STEP_FLAG, step + "", "--" + WORDSIZE_FLAG);
    }

    final CountParams countParams = new CountParams(output, 1/*topn*/, 1/*min*/, false);

    final BuildSearchParamsBuilder builder = BuildSearchParams.builder()
        .mode(pm).count(countParams).uniqueWords(mFlags.isSet(UNIQUE_WORDS));

    if (mFlags.isSet(INPUT_FLAG)) {
      final File subject = (File) mFlags.getValue(INPUT_FLAG);
      final SequenceParams subjectParams = SequenceParams.builder().directory(subject).mode(pm.subjectMode()).create();
      final int idBits = MathUtils.ceilPowerOf2Bits(subjectParams.numberSequences());
      Diagnostic.developerLog("idBits=" + idBits + " subjectParams.numberSequences=" + subjectParams.numberSequences());
      final long totalLength = subjectParams.reader().totalLength();
      mUsageMetric.setMetric(totalLength);
      final long size = BuildParams.size(totalLength, subjectParams.numberSequences(), nwindow, nstep, 1, pm.subjectMode().codeIncrement());
      final BuildParams buildParams = BuildParams.builder()
          .windowSize(nwindow)
          .stepSize(nstep)
          .idBits(idBits)
          .compressHashes(true)
          .sequences(subjectParams)
          .createBitVector(false)
          .spaceEfficientButUnsafe(true)
          .ideal(true)
          .size(size)
          .create();
      builder.build(buildParams);
    } else {
      final File listFile = (File) mFlags.getValue(CommonFlags.INPUT_LIST_FLAG);
      int numberErrorLines = 0;
      try (BufferedReader inputList = new BufferedReader(new FileReader(listFile))) {
        final LongRange readerRestriction = mFlags.isSet(MAX_READS_FLAG) ? new LongRange(0, (Integer) mFlags.getValue(MAX_READS_FLAG)) : LongRange.NONE;
        final Map<String, List<SequenceParams>> labelMap = new HashMap<>();
        final ArrayList<String> orderedLabels = new ArrayList<>();
        long length = 0;
        long numSeqs = 0;
        String line;
        while ((line = inputList.readLine()) != null) {
          final int spaceIndex = line.indexOf(' ');
          if (spaceIndex == -1 || spaceIndex == 0 || spaceIndex == line.length() - 1) {
            continue;
          }
          final String label = line.substring(0, spaceIndex).trim();
          final String fileName = line.substring(spaceIndex).trim();
          final File currentFile = new File(fileName);
          if (!currentFile.exists()) {
            if (numberErrorLines < 10) {
              Diagnostic.error(ErrorType.INFO_ERROR, "The specified SDF, \"" + currentFile.getPath() + "\", does not exist.");
            }
            numberErrorLines++;
            continue;
          } else if (!currentFile.isDirectory()) {
            if (numberErrorLines < 10) {
              Diagnostic.error(ErrorType.INFO_ERROR, "The specified file, \"" + currentFile.getPath() + "\", is not an SDF.");
            }
            numberErrorLines++;
            continue;
          }
          final List<SequenceParams> paramsList;
          if (labelMap.containsKey(label)) {
            paramsList = labelMap.get(label);
          } else {
            paramsList = new ArrayList<>();
            labelMap.put(label, paramsList);
            orderedLabels.add(label);
          }

          if (ReaderUtils.isPairedEndDirectory(currentFile)) {
            final LongRange resolvedRange = resolveRange(currentFile, ReaderUtils.getLeftEnd(currentFile), readerRestriction);
            if (resolvedRange == null) {
              continue;
            }
            final SequenceParams subjectLeftParams = SequenceParams.builder().readerRestriction(resolvedRange).directory(ReaderUtils.getLeftEnd(currentFile)).mode(pm.subjectMode()).create();
            final SequenceParams subjectRightParams = SequenceParams.builder().readerRestriction(resolvedRange).directory(ReaderUtils.getRightEnd(currentFile)).mode(pm.subjectMode()).create();
            paramsList.add(subjectLeftParams);
            length += subjectLeftParams.reader().totalLength();
            numSeqs += subjectLeftParams.reader().numberSequences();
            paramsList.add(subjectRightParams);
            length += subjectRightParams.reader().totalLength();
            numSeqs += subjectRightParams.reader().numberSequences();
          } else {
            final LongRange resolvedRange = resolveRange(currentFile, currentFile, readerRestriction);
            if (resolvedRange == null) {
              continue;
            }
            final SequenceParams subjectParams = SequenceParams.builder().readerRestriction(resolvedRange).directory(currentFile).mode(pm.subjectMode()).create();
            paramsList.add(subjectParams);
            length += subjectParams.reader().totalLength();
            numSeqs += subjectParams.reader().numberSequences();
          }
        }
        if (numberErrorLines > 0) {
          throw new InvalidParamsException(ErrorType.INFO_ERROR, "There were " + numberErrorLines + " incorrect input list lines.");
        }
        if (labelMap.isEmpty()) {
          throw new InvalidParamsException(ErrorType.INFO_ERROR, "The input list file contained no target files.");
        }
        final List<Pair<String, List<SequenceParams>>> subjects = new ArrayList<>();
        assert orderedLabels.size() == labelMap.size();
        for (final String label : orderedLabels) {
          subjects.add(new Pair<>(label, labelMap.get(label)));
        }
        builder.sequences(subjects);
        final SequenceParams dummySubjectParams = SequenceParams.builder().region(new HashingRegion(0, subjects.size())).mode(pm.subjectMode()).create();
        final long size = BuildParams.size(length, numSeqs, nwindow, nstep, 1, pm.subjectMode().codeIncrement());
        final int idBits = MathUtils.ceilPowerOf2Bits(orderedLabels.size());
        Diagnostic.developerLog("idBits=" + idBits + " orderedLabels.size=" + orderedLabels.size());
        final BuildParams buildParams = BuildParams.builder()
          .windowSize(nwindow)
          .stepSize(nstep)
          .size(size)
          .idBits(idBits)
          .compressHashes(true)
          .sequences(dummySubjectParams)
          .createBitVector(false)
          .spaceEfficientButUnsafe(true)
          .ideal(true)
          .create();
        builder.build(buildParams);
      }
    }

    return builder.create();
  }

  protected static LongRange resolveRange(File currentFile, File sdfFile, LongRange region) throws IOException {
    final IndexFile index = new IndexFile(sdfFile);
    if (index.getNumberSequences() == 0) {
      Diagnostic.warning("The SDF \"" + currentFile.getPath() + "\" contains no sequences");
      return null;
    }
    return SequencesReaderFactory.resolveRange(index, region);
  }

  @Override
  protected IORunnable task(final BuildSearchParams params, final OutputStream out) {
    return new SimilarityTask(params, out, mUsageMetric);
  }

}

