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
package com.rtg.protein;

import static com.rtg.launcher.BuildCommon.RESOURCE;
import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.launcher.CommonFlags.READS_FLAG;
import static com.rtg.launcher.CommonFlags.TEMPLATE_FLAG;
import static com.rtg.launcher.CommonFlags.TEMP_DIR;
import static com.rtg.ngs.MapFlags.COMPRESS_HASHES_FLAG;
import static com.rtg.ngs.MapFlags.DEFAULT_TOP_N;
import static com.rtg.ngs.MapFlags.TEMP_FILES_COMPRESSED;
import static com.rtg.ngs.MapFlags.WORDSIZE_FLAG;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.index.params.CreateParams;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.ParamsCli;
import com.rtg.launcher.ParamsTask;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.mode.SequenceMode;
import com.rtg.ngs.MapFlags;
import com.rtg.ngs.MapParamsHelper;
import com.rtg.ngs.NameParams;
import com.rtg.ngs.NgsFilterParams;
import com.rtg.ngs.NgsMaskParams;
import com.rtg.ngs.NgsMaskParamsProtein;
import com.rtg.ngs.NgsOutputParams;
import com.rtg.ngs.NgsOutputParamsBuilder;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.ngs.OutputFilter;
import com.rtg.ngs.SamSequenceReaderParams;
import com.rtg.reader.FormatCli;
import com.rtg.reader.ReaderUtils;
import com.rtg.report.MapSummaryReport;
import com.rtg.report.MapXSummaryReport;
import com.rtg.report.ReportType;
import com.rtg.usage.UsageMetric;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.MathUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.Timer;

/**
 * Module for doing protein matching
 *
 */
public class MapXCli extends ParamsCli<NgsParams> {

  private static final String UNFILTERED_FLAG = "all-hits";

  private static final String GAPS_FLAG = "gaps";
  static final String GAP_LENGTH_FLAG = "gap-length";
  private static final String MISMATCHES_FLAG = "mismatches";
  static final String MATRIX_FLAG = "matrix";
  private static final String READ_CACHE_FLAG = "Xenable-read-cache";
  static final String PRE_FILTER_ALGORITHM = "Xprefilter-algorithm";
  private static final String PRE_FILTER_MIN_SCORE = "Xprefilter-min-score";
  private static final String PRE_FILTER_MIN_OVERLAP = "Xprefilter-min-overlap";
  static final String MIN_IDENTITY_FLAG = "min-identity";
  static final String MAX_ESCORE_FLAG = "max-e-score";
  static final String MIN_BITSCORE_FLAG = "min-bit-score";
  static final String MIN_READ_LENGTH = "min-dna-read-length";
  private static final String OUTPUT_READ_NAMES_FLAG = "read-names";
  private static final String SUPPRESS_PROTEIN_OUTPUT_FLAG = "suppress-protein";
  private static final String XMETA_CHUNK_LENGTH = "Xmeta-chunk-length";
  private static final String XMETA_CHUNK_OVERLAP = "Xmeta_chunk-overlap";

  private static final String XDONT_MERGE_ALIGNMENT_RESULTS = "Xdont-merge-alignment-result";

  static final HashMap<String, OutputFilter> FILTERS = new HashMap<>();
  private static final ArrayList<String> MATRICES = new ArrayList<>();

  static {
    FILTERS.put(MapFlags.TOPEQUAL, OutputFilter.PROTEIN_TOPEQUAL);
    FILTERS.put(MapFlags.TOPN, OutputFilter.PROTEIN_TOPN);
    MATRICES.add("blosum45");
    MATRICES.add("blosum62");
    MATRICES.add("blosum80");
  }

  /** Max alignment score */
  public static final String MAX_ALIGNMENT_SCORE = "max-alignment-score";

  @Override
  public String moduleName() {
    return "mapx";
  }

  @Override
  public String description() {
    return "translated protein search";
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(OUTPUT_FLAG);
  }

  @TestClass(value = "com.rtg.protein.MapXValidatorTest")
  private static class MapXValidator implements Validator {

    @Override
    public boolean isValid(final CFlags flags) {
      final String format = flags.isSet(FormatCli.FORMAT_FLAG) ? flags.getValue(FormatCli.FORMAT_FLAG).toString().toLowerCase(Locale.getDefault()) : FormatCli.SDF_FORMAT;
      final boolean sdf = format.equals(FormatCli.SDF_FORMAT);
      if (!CommonFlags.validateOutputDirectory(flags) || !CommonFlags.validateReads(flags, sdf)) {
        return false;
      }
      if (!CommonFlags.validateTemplate(flags) || !CommonFlags.validateThreads(flags) || !MapFlags.checkPercentRepeatFrequency(flags)) {
        return false;
      }

      if (sdf && ReaderUtils.isPairedEndDirectory((File) flags.getValue(READS_FLAG))) {
        flags.setParseMessage("Paired end data not supported");
        return false;
      }
      final int a = (Integer) flags.getValue(MISMATCHES_FLAG);
      final int b = (Integer) flags.getValue(GAPS_FLAG);
      final int c = (Integer) flags.getValue(GAP_LENGTH_FLAG);
      if (flags.isSet(WORDSIZE_FLAG)) {
        if (!flags.checkInRange(WORDSIZE_FLAG, 1, 12)) {
          return false;
        }
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
      if (!flags.checkInRange(MIN_IDENTITY_FLAG, 0, 100)) {
        return false;
      }
      if (!flags.checkNand(MAX_ESCORE_FLAG, MIN_BITSCORE_FLAG)) {
        return false;
      }
      if (!flags.checkInRange(MAX_ESCORE_FLAG, 0, Double.MAX_VALUE)) {
        return false;
      }
      if (!flags.checkInRange(MIN_BITSCORE_FLAG, 0, Double.MAX_VALUE)) {
        return false;
      }
      if (!flags.checkInRange(PRE_FILTER_MIN_SCORE, 0, 100)) {
        return false;
      }
      if (!flags.checkInRange(PRE_FILTER_MIN_OVERLAP, 0, 100)) {
        return false;
      }
      if (!flags.checkInRange(PRE_FILTER_ALGORITHM, -10, 10)) {
        return false;
      }
      final int metachunklength;
      if (flags.isSet(XMETA_CHUNK_LENGTH)) {
        metachunklength = (Integer) flags.getValue(XMETA_CHUNK_LENGTH);
      } else {
        metachunklength = ProteinReadIndexer.DEFAULT_META_CHUNK_LENGTH;
      }
      final int metachunkoverlap;
      if (flags.isSet(XMETA_CHUNK_OVERLAP)) {
        metachunkoverlap = (Integer) flags.getValue(XMETA_CHUNK_OVERLAP);
      } else {
        metachunkoverlap = metachunklength / 2;
      }
      if (metachunklength < 1 || metachunklength > 63) {
        flags.setParseMessage("--" + XMETA_CHUNK_LENGTH + " must be set between " + 1 + " and " + 63);
        return false;
      }
      if (metachunkoverlap >= metachunklength || metachunkoverlap < 0) {
        flags.setParseMessage("--" + XMETA_CHUNK_OVERLAP + " must be positive and less than " + XMETA_CHUNK_LENGTH);
        return false;
      }
      if (!flags.checkInRange(MapFlags.MAX_TOP_RESULTS_FLAG, 1, 250)) {
        return false;
      }
      return true;
    }
  }

  @Override
  protected NgsParams makeParams() throws InvalidParamsException, IOException {
    final NgsParamsBuilder ngsParamsBuilder = NgsParams.builder();
    final NameParams nameParams = new NameParams(mFlags.isSet(OUTPUT_READ_NAMES_FLAG), false);
    MapParamsHelper.initReaders(ngsParamsBuilder, mFlags, nameParams, false, SequenceMode.PROTEIN, SequenceMode.TRANSLATED, new SamSequenceReaderParams(false, false));
    final ISequenceParams reads = ngsParamsBuilder.buildFirstParams();
    final ISequenceParams templ = ngsParamsBuilder.searchParams();
    final int readLen = (int) reads.reader().maxLength();
    final int minReadLen = (int) reads.reader().minLength();
    final NgsMaskParams mask = createMaskParams(readLen);
    ngsParamsBuilder.maskParams(mask);
    final long mapXMinLength;
    if (mFlags.isSet(MIN_READ_LENGTH)) {
      mapXMinLength = (Long) mFlags.getValue(MIN_READ_LENGTH);
    } else {
      mapXMinLength = 3 * (mask.getWordSize() + mask.getSubstitutions() + 1);
    }
    ngsParamsBuilder.mapXMinLength(mapXMinLength);
    if (readLen < mapXMinLength) {
      throw new NoTalkbackSlimException("All reads are shorter than the minimum DNA read length " + mapXMinLength);
    }
    if (minReadLen < mapXMinLength) {
      Diagnostic.warning("The read set contains reads which are shorter than the minimum DNA read length " + mapXMinLength  + " which will be ignored");
    }

    final NgsFilterParams.NgsFilterParamsBuilder filterBuild = NgsFilterParams.builder();
    buildFilterParams(mFlags, filterBuild);

    if (mFlags.isSet(MIN_BITSCORE_FLAG)) {
      filterBuild.minBitScore((Double) mFlags.getValue(MIN_BITSCORE_FLAG));
      filterBuild.maxEScore(Double.MAX_VALUE);
    } else {
      filterBuild.maxEScore((Double) mFlags.getValue(MAX_ESCORE_FLAG));
    }
    final NgsFilterParams filter = filterBuild.create();

    final NgsOutputParamsBuilder outBuild = NgsOutputParams.builder();
    outBuild.outputDir((File) mFlags.getValue(OUTPUT_FLAG)).filterParams(filter).sorted(mFlags.isSet(CommonFlags.SORT_FLAG));
    outBuild.tempFilesDir((File) mFlags.getValue(TEMP_DIR));
    outBuild.outputUnmapped(!mFlags.isSet(MapFlags.NO_UNMAPPED));
    if (mFlags.isSet(XDONT_MERGE_ALIGNMENT_RESULTS)) {
      outBuild.mergeAlignmentResults(false);
    }
    outBuild.mergeMatchResults(false);
    outBuild.outputReadNames(mFlags.isSet(OUTPUT_READ_NAMES_FLAG));
    outBuild.outputProteinSequences(!mFlags.isSet(SUPPRESS_PROTEIN_OUTPUT_FLAG));
    ngsParamsBuilder.compressHashes((Boolean) mFlags.getValue(COMPRESS_HASHES_FLAG));
    ngsParamsBuilder.outputParams(outBuild.create());
    ngsParamsBuilder.useLongReadMapping(false);
    /*
    if ((readLen / 3 - 1) > 64) {
      throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "Input SDF contains reads that will translate to greater than 64 amino acids");
      //ngsParamsBuilder.useLongReadMapping(true);
    }*/
    ngsParamsBuilder.buildFirstParams(reads);
    ngsParamsBuilder.searchParams(templ);
    if (mFlags.getFlag(MATRIX_FLAG) != null) {
      ngsParamsBuilder.proteinScoringMatrix(new ProteinScoringMatrix((String) mFlags.getValue(MATRIX_FLAG)));
    }
    ngsParamsBuilder.enableProteinReadCache(mFlags.isSet(READ_CACHE_FLAG));

    MapParamsHelper.populateProportionalRepeat(mFlags, ngsParamsBuilder);

    ngsParamsBuilder.numberThreads(CommonFlags.parseThreads((Integer) mFlags.getValue(CommonFlags.THREADS_FLAG)));
    final int metaChunkSize;
    if (mFlags.isSet(XMETA_CHUNK_LENGTH)) {
      metaChunkSize = (Integer) mFlags.getValue(XMETA_CHUNK_LENGTH);
    } else {
      metaChunkSize = ProteinReadIndexer.DEFAULT_META_CHUNK_LENGTH;
    }
    final int metaChunkOverlap;
    if (mFlags.isSet(XMETA_CHUNK_OVERLAP)) {
      metaChunkOverlap = (Integer) mFlags.getValue(XMETA_CHUNK_OVERLAP);
    } else {
      metaChunkOverlap = metaChunkSize / 2;
    }
    ngsParamsBuilder.mapXMetaChunkSize(metaChunkSize);
    ngsParamsBuilder.mapXMetaChunkOverlap(metaChunkOverlap);
    final NgsParams localParams = ngsParamsBuilder.create();
    if (localParams.outputParams().outputReadNames() && !localParams.buildFirstParams().reader().hasNames()) {
      throw new InvalidParamsException("Names not present in SDF and read names requested");
    }
    localParams.globalIntegrity();
    return localParams;
  }

  /**
   * Build filter parameters
   * @param flags command line flags
   * @param builder filter builder
   */
  static void buildFilterParams(final CFlags flags, final NgsFilterParams.NgsFilterParamsBuilder builder) {
    final OutputFilter outputFilter = flags.isSet(UNFILTERED_FLAG) ? OutputFilter.PROTEIN_ALL_HITS : FILTERS.get(((String) flags.getValue(CommonFlags.OUTPUT_FILTER_FLAG)).toLowerCase(Locale.getDefault()));

    final int maxTopResults;
    if (flags.isSet(MapFlags.MAX_TOP_RESULTS_FLAG)) {
      maxTopResults = (Integer) flags.getValue(MapFlags.MAX_TOP_RESULTS_FLAG);
    } else {
      maxTopResults = DEFAULT_TOP_N;
    }
    final int topN;
    if (flags.isSet(MapFlags.TOPN_RESULTS_FLAG)) {
      topN = (Integer) flags.getValue(MapFlags.TOPN_RESULTS_FLAG);
    } else {
      // mapx doesn't use TOPN_FLAG, only MAX_TOP_RESULTS_FLAG
      if (flags.getFlag(MapFlags.TOPN_RESULTS_FLAG) == null && flags.isSet(MapFlags.MAX_TOP_RESULTS_FLAG)) {
        topN = (Integer) flags.getValue(MapFlags.MAX_TOP_RESULTS_FLAG);
      } else {
        topN = Math.max(maxTopResults, DEFAULT_TOP_N);
      }
    }

    final Integer minIdentity = (Integer) flags.getValue(MIN_IDENTITY_FLAG);
    final Integer minScore = flags.isSet(PRE_FILTER_MIN_SCORE) ? (Integer) flags.getValue(PRE_FILTER_MIN_SCORE) : Integer
        .valueOf((int) (0.8 * minIdentity));
    final Integer minOverlap = flags.isSet(PRE_FILTER_MIN_OVERLAP) ? (Integer) flags.getValue(PRE_FILTER_MIN_OVERLAP) : minIdentity;

    final boolean zip = !flags.isSet(CommonFlags.NO_GZIP);
    builder.outputFilter(outputFilter)
    .zip(zip)
    .topN(topN)
    .maxTopResults(maxTopResults)
    .exclude(flags.isSet(CommonFlags.EXCLUDE_FLAG))
    .useids(flags.isSet(CommonFlags.USEIDS_FLAG))
    .errorLimit(flags.isSet(MapFlags.XSCORE_INDEL) ? (Integer) flags.getValue(MapFlags.XSCORE_INDEL) : MapFlags.MAX_SCORE)
    .matedMaxMismatches((IntegerOrPercentage) flags.getValue(MapXCli.MAX_ALIGNMENT_SCORE))
    .minIdentity(minIdentity)
    .preFilterAlgorithm((Integer) flags.getValue(PRE_FILTER_ALGORITHM))
    .preFilterMinScore(minScore)
    .preFilterMinOverlap(minOverlap);
  }

  private NgsMaskParams createMaskParams(final int readLength) {
    final int w = (Integer) mFlags.getValue(WORDSIZE_FLAG);
    final int subs = (Integer) mFlags.getValue(MISMATCHES_FLAG);
    final int i = (Integer) mFlags.getValue(GAPS_FLAG);
    final int l = (Integer) mFlags.getValue(GAP_LENGTH_FLAG);
    final int s = Math.max(subs, i);
    MapFlags.validateMaskParams(readLength / 3 - 1, w, s, i, l);
    return new NgsMaskParamsProtein(w, s, i, l);
  }

  @Override
  protected ParamsTask<?, ?> task(final NgsParams params, final OutputStream out) {
    final UsageMetric usageMetric = mUsageMetric == null ? new UsageMetric() : mUsageMetric; //create when null to cover some testing
    return new MapXTask(params, out, usageMetric);
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  private void initFlags(final CFlags flags) {
    flags.registerExtendedHelp();
    flags.setDescription("Searches translated read data sets of defined length (e.g. 100 bp reads) against protein databases or translated nucleotide sequences.");
    CommonFlagCategories.setCategories(flags);
    flags.registerRequired('i', READS_FLAG, File.class, "SDF|FILE", "query read sequences").setCategory(INPUT_OUTPUT);
    flags.registerRequired('t', TEMPLATE_FLAG, File.class, "SDF", "SDF containing database to search").setCategory(INPUT_OUTPUT);
    CommonFlags.initOutputDirFlag(flags);

    // No Paired End input for MapX
    final Flag formatFlag = flags.registerOptional('F', FormatCli.FORMAT_FLAG, String.class, "FORMAT", "input format for reads", FormatCli.SDF_FORMAT).setCategory(INPUT_OUTPUT);
    formatFlag.setParameterRange(new String[] {FormatCli.SDF_FORMAT, FormatCli.FASTA_FORMAT, FormatCli.FASTQ_FORMAT, FormatCli.SAM_SE_FORMAT});

    final Flag filter = flags.registerOptional('f', CommonFlags.OUTPUT_FILTER_FLAG, String.class, "STRING", "output filter", "topn");
    filter.setCategory(REPORTING);
    filter.setParameterRange(FILTERS.keySet());
    flags.registerOptional(XDONT_MERGE_ALIGNMENT_RESULTS, "does not concat alignment files").setCategory(UTILITY);
    flags.registerOptional('n', MapFlags.MAX_TOP_RESULTS_FLAG, Integer.class, "int", "maximum number of topn/topequals results output per read",
      DEFAULT_TOP_N).setCategory(REPORTING);
    final Flag matrix = flags.registerOptional(MATRIX_FLAG, String.class, "string", "name of the scoring matrix used during alignment", "blosum62");
    matrix.setCategory(SENSITIVITY_TUNING);
    matrix.setParameterRange(MATRICES);
    flags.registerOptional(PRE_FILTER_ALGORITHM, Integer.class, "int", "pre-filter algorithm", -3).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(PRE_FILTER_MIN_SCORE, Integer.class, "int", "pre-filter minimum score").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(PRE_FILTER_MIN_OVERLAP, Integer.class, "int", "pre-filter minimum overlap").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('P', MIN_IDENTITY_FLAG, Integer.class, "int", "minimum percent identity at output", 60).setCategory(REPORTING);
    flags.registerOptional('E', MAX_ESCORE_FLAG, Double.class, "float", "maximum e-score at output", 10.0).setCategory(REPORTING);
    flags.registerOptional('B', MIN_BITSCORE_FLAG, Double.class, "float", "minimum bit score at output").setCategory(REPORTING);
    flags.registerOptional(COMPRESS_HASHES_FLAG, Boolean.class, "BOOL", "compress hashes in indexes", true).setCategory(UTILITY);
    flags.registerOptional(TEMP_DIR, File.class, "DIR", "directory used for temporary files (Defaults to output directory)").setCategory(UTILITY);
    flags.registerOptional(TEMP_FILES_COMPRESSED, Boolean.class, "BOOL", "gzip temporary files", true).setCategory(UTILITY);
    flags.registerOptional(MapFlags.NO_UNMAPPED, "do not output unmapped reads").setCategory(UTILITY);
    flags.registerOptional('w', WORDSIZE_FLAG, Integer.class, "int", "word size", 7).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('a', MISMATCHES_FLAG, Integer.class, "INT", "guaranteed minimum number of identical mismatches which will be detected", 1).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('b', GAPS_FLAG, Integer.class, "INT", "guaranteed minimum number of gaps which will be detected (if this is larger than the minimum number of mismatches then the minimum number of mismatches is increased to the same value)", 0).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('c', GAP_LENGTH_FLAG, Integer.class, "INT", "guaranteed number of positions that will be detected in a single gap", 1).setCategory(SENSITIVITY_TUNING);
    CommonFlags.initThreadsFlag(mFlags);
    CommonFlags.initNoGzip(flags);
    flags.registerOptional('e', MAX_ALIGNMENT_SCORE, IntegerOrPercentage.class, "INT", "maximum alignment score at output (as absolute value or percentage of read length in protein space)", IntegerOrPercentage.valueOf("30%")).setCategory(REPORTING);
    flags.registerOptional(MapFlags.REPEAT_FREQUENCY_FLAG, IntegerOrPercentage.class, "INT", RESOURCE.getString("REPEAT_FREQUENCY_DESC"), IntegerOrPercentage.valueOf("95%")).setCategory(SENSITIVITY_TUNING);

    flags.registerOptional(OUTPUT_READ_NAMES_FLAG, "use read name in output instead of read id (Uses more RAM)").setCategory(UTILITY);
    flags.registerOptional(SUPPRESS_PROTEIN_OUTPUT_FLAG, "suppress output of sequence protein information").setCategory(UTILITY);

    flags.registerOptional(READ_CACHE_FLAG, "enable protein read cache").setCategory(UTILITY);

    flags.registerOptional(UNFILTERED_FLAG, "output all alignments meeting thresholds instead of applying topn/topequals N limits").setCategory(REPORTING);
    flags.registerOptional(MIN_READ_LENGTH, Long.class, "INT", "minimum read length in nucleotides. Shorter reads will be ignored. (Defaults to 3 * (w + a + 1))").setCategory(SENSITIVITY_TUNING);

    flags.registerOptional(XMETA_CHUNK_LENGTH, Integer.class, "INT", "how large to make long read meta chunks. (Defaults to " + ProteinReadIndexer.DEFAULT_META_CHUNK_LENGTH + ")").setCategory(SENSITIVITY_TUNING);
    flags.registerOptional(XMETA_CHUNK_OVERLAP, Integer.class, "INT", "how much overlap to have in long read meta chunks. (Defaults to meta chunk length / 2)").setCategory(SENSITIVITY_TUNING);

    CommonFlags.initReadRange(mFlags);

    flags.setValidator(new MapXValidator());
  }

  @TestClass("com.rtg.protein.MapXFunctionalTest")
  private static class MapXTask extends ParamsTask<NgsParams, MapXStatistics> {

    MapXTask(final NgsParams params, final OutputStream reportStream, final UsageMetric usageMetric) {
      super(params, reportStream, new MapXStatistics(params.directory()), usageMetric);
    }

    private void execNgs(final NgsParams params) throws IOException {
      //make all the components we need
      final long indexSize = params.buildFirstParams().numberSequences() * params.buildFirstParams().mode().numberFrames();
      assert params.searchParams().numberSequences() < Integer.MAX_VALUE : indexSize;
      final OutputFilter filter = params.outputParams().outFilter();
      try (OutputProcessor outProcessor = filter.makeProcessor(params, mStatistics)) {
        final int numChunks = ProteinReadIndexer.countMetaChunks((int) params.buildFirstParams().maxLength(), params.mapXMetaChunkSize(), params.mapXMetaChunkOverlap());
        final long numValues = params.buildFirstParams().numberSequences() * params.buildFirstParams().mode().numberFrames() * numChunks;
        assert numValues < Integer.MAX_VALUE : numValues;
        final int bitValues = MathUtils.ceilPowerOf2Bits(numValues);

        final SequenceLengthBuckets buckets = new SequenceLengthBuckets(params.buildFirstParams().reader(), params.mapXMinReadLength());
        final CreateParams.CreateParamsBuilder builder = new CreateParams.CreateParamsBuilder();
        builder.valueBits(bitValues)
                .compressHashes(params.compressHashes());

        mUsageMetric.setMetric(ProteinReadIndexer.indexThenSearchProteinReads(params, outProcessor, builder, buckets, numValues));

        outProcessor.finish();
      }
      final MapSummaryReport mr = new MapXSummaryReport(); //this filter params ought to be ignored, since the data exists...
      mr.setCommandLines(Collections.singletonList(CommandLine.getCommandLine()));
      mr.setParams(mParams);
      final List<File> samples = new ArrayList<>();
      final File left = mParams.buildFirstParams().reader().path();
      // Don't have to worry about paired end in MAPX
      assert mParams.buildSecondParams() == null;
      samples.add(left);
      mr.setSampleFiles(samples);
      mr.setTitle("Protein mapping report");
      mr.generateReport(ReportType.HTML, mParams.directory(), mParams.directory());
    }

    @Override
    protected void exec() throws IOException {
      try (NgsParams proteinParams = mParams) {
        final Timer firstTimer = new Timer("Matching_Phase");

        firstTimer.start();
        execNgs(proteinParams);
        firstTimer.stop();
        firstTimer.log();
      } finally {
        mParams.close();
      }
    }
  }

}
