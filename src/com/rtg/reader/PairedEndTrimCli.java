/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.reader;

import static com.rtg.launcher.CommonFlags.MIN_READ_LENGTH;
import static com.rtg.launcher.CommonFlags.NO_GZIP;
import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.launcher.CommonFlags.QUALITY_FLAG;
import static com.rtg.reader.FastqTrim.BATCH_SIZE;
import static com.rtg.reader.Sdf2Fasta.INTERLEAVE;
import static com.rtg.sam.SamFilterOptions.SUBSAMPLE_FLAG;
import static com.rtg.util.cli.CommonFlagCategories.FILTERING;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.Collections;
import java.util.concurrent.FutureTask;
import java.util.function.Function;
import java.util.function.Predicate;

import com.rtg.alignment.SingleIndelSeededEditDistance;
import com.rtg.alignment.UnidirectionalAdaptor;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.ngs.MapFlags;
import com.rtg.ngs.MapParamsHelper;
import com.rtg.ngs.NgsParams;
import com.rtg.sam.SamFilterOptions;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.Timer;
import com.rtg.util.io.BaseFile;
import com.rtg.util.io.FileUtils;

/**
 */
public class PairedEndTrimCli extends AbstractCli {

  private static final String PETRIM_MODULE_NAME = "petrim";


  private static final String RIGHT = "right";
  private static final String LEFT = "left";

  private static final String MIDPOINT_TRIM = "midpoint-trim";

  private static final String MIN_IDENTITY = "min-identity";
  private static final String MIN_OVERLAP = "min-overlap-length";

  private static final String LEFT_PROBE_LENGTH = "left-probe-length";
  private static final String RIGHT_PROBE_LENGTH = "right-probe-length";

  private static final String DISCARD_EMPTY_READS = "discard-empty-reads";
  private static final String DISCARD_EMPTY_PAIRS = "discard-empty-pairs";

  private static final String VERBOSE = "Xverbose";


  @Override
  protected void initFlags() {
    mFlags.setDescription(StringUtils.sentencify(description()));
    initFlags(mFlags);
  }

  static void initFlags(CFlags flags) {
    flags.registerExtendedHelp();
    CommonFlagCategories.setCategories(flags);
    flags.registerRequired('l', LEFT, File.class, CommonFlags.FILE, "left input FASTQ file (AKA R1)").setCategory(INPUT_OUTPUT);
    flags.registerRequired('r', RIGHT, File.class, CommonFlags.FILE, "right input FASTQ file (AKA R2)").setCategory(INPUT_OUTPUT);
    flags.registerRequired('o', OUTPUT_FLAG, File.class, CommonFlags.FILE, "output filename prefix. Use '-' to write to standard output").setCategory(INPUT_OUTPUT);
    CommonFlags.initQualityFormatFlag(flags);
    MapFlags.initAlignerPenaltyFlags(flags);
    CommonFlags.initMinReadLength(flags);
    flags.registerOptional('L', MIN_OVERLAP, Integer.class, CommonFlags.INT, "minimum number of bases in overlap to trigger overlap trimming", 25).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('P', MIN_IDENTITY, Integer.class, CommonFlags.INT, "minimum percent identity in overlap to trigger overlap trimming", 90).setCategory(SENSITIVITY_TUNING);

    flags.registerOptional('m', MIDPOINT_TRIM, "if set, trim both reads to midpoint of overlap region").setCategory(FILTERING);
    flags.registerOptional(LEFT_PROBE_LENGTH, Integer.class, CommonFlags.INT, "assume R1 starts with probes this long, and trim R2 bases that overlap into this", 0).setCategory(FILTERING);
    flags.registerOptional(RIGHT_PROBE_LENGTH, Integer.class, CommonFlags.INT, "assume R2 starts with probes this long, and trim R1 bases that overlap into this", 0).setCategory(FILTERING);
    flags.registerOptional(DISCARD_EMPTY_PAIRS, "if set, discard pairs where both reads have zero length (after any trimming)").setCategory(FILTERING);
    flags.registerOptional(DISCARD_EMPTY_READS, "if set, discard pairs where either read has zero length (after any trimming)").setCategory(FILTERING);

    flags.registerOptional(INTERLEAVE, "interleave paired data into a single output file. Default is to split to separate output files").setCategory(UTILITY);
    CommonFlags.initThreadsFlag(flags);
    SamFilterOptions.registerSubsampleFlags(flags);
    flags.registerOptional(BATCH_SIZE, Integer.class, CommonFlags.INT, "number of pairs to process per batch", 10000).setCategory(UTILITY);
    flags.registerOptional(VERBOSE, "dump read alignment information to stderr").setCategory(UTILITY);
    CommonFlags.initNoGzip(flags);

    flags.setValidator(new FlagsValidator());
  }

  static class FlagsValidator implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      final File baseOutput = (File) flags.getValue(OUTPUT_FLAG);
      final boolean gzip = !flags.isSet(NO_GZIP);
      final BaseFile baseFile = FastqUtils.baseFile(baseOutput, gzip);
      if (flags.isSet(INTERLEAVE)) {
        if (!CommonFlags.validateOutputFile(flags, baseFile.suffixedFile(""))) {
          return false;
        }

      } else {
        if (FileUtils.isStdio(baseOutput)) {
          flags.setParseMessage("Sending non-interleaved paired-end data to stdout is not supported.");
          return false;
        }
        if (!(CommonFlags.validateOutputFile(flags, baseFile.suffixedFile("_1")) && CommonFlags.validateOutputFile(flags, baseFile.suffixedFile("_2")))) {
          return false;
        }
      }
      return CommonFlags.validateInputFile(flags, LEFT)
        && CommonFlags.validateInputFile(flags, RIGHT)
        && flags.checkInRange(BATCH_SIZE, 1, Integer.MAX_VALUE)
        && flags.checkInRange(MIN_OVERLAP, 1, Integer.MAX_VALUE)
        && flags.checkInRange(MIN_IDENTITY, 1, 100)
        && flags.checkInRange(SUBSAMPLE_FLAG, 0.0, 1.0)
        && flags.checkInRange(MIN_READ_LENGTH, 0, Integer.MAX_VALUE)
        && flags.checkInRange(LEFT_PROBE_LENGTH, 0, Integer.MAX_VALUE)
        && flags.checkInRange(RIGHT_PROBE_LENGTH, 0, Integer.MAX_VALUE)
        && flags.checkNand(DISCARD_EMPTY_PAIRS, DISCARD_EMPTY_READS);
    }
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final boolean gzip = !mFlags.isSet(NO_GZIP);
    final File baseOutput = (File) mFlags.getValue(OUTPUT_FLAG);
    final boolean interleavePaired = mFlags.isSet(INTERLEAVE);
    final int batchSize = (Integer) mFlags.getValue(BATCH_SIZE);
    final PairAlignmentStats stats = new PairAlignmentStats();
    final QualityFormat encoding = FastqTrim.qualityFlagToFastQScoreType((String) mFlags.getValue(QUALITY_FLAG));
    // All trimming and aligning is done in separate threads from reading
    final int threads = mFlags.isSet(VERBOSE) ? 1 : CommonFlags.parseThreads((Integer) mFlags.getValue(CommonFlags.THREADS_FLAG));
    try (final SequenceDataSource r1fq = new FastqSequenceDataSource(Collections.singletonList(FileUtils.createInputStream((File) mFlags.getValue(LEFT), true)), encoding);
         final SequenceDataSource r2fq = new FastqSequenceDataSource(Collections.singletonList(FileUtils.createInputStream((File) mFlags.getValue(RIGHT), true)), encoding)) {

      final Timer t = new Timer("FastqPairTrimmer");
      t.start();

      final BaseFile baseFile = FastqUtils.baseFile(baseOutput, gzip);
      final FastqWriter left;
      final FastqWriter right;
      if (interleavePaired) {
        left = new FastqWriter(new OutputStreamWriter(FileUtils.createOutputStream(baseFile, "")));
        right = left;
      } else {
        if (FileUtils.isStdio(baseOutput)) {
          throw new NoTalkbackSlimException("Sending non-interleaved paired-end data to stdout is not supported.");
        }
        left = new FastqWriter(new OutputStreamWriter(FileUtils.createOutputStream(baseFile, "_1")));
        right = new FastqWriter(new OutputStreamWriter(FileUtils.createOutputStream(baseFile, "_2")));
      }
      final Predicate<FastqPair> p;
      if (mFlags.isSet(DISCARD_EMPTY_READS)) {
        p = o -> o.r1().length() != 0 && o.r2().length() != 0;
      } else if (mFlags.isSet(DISCARD_EMPTY_PAIRS)) {
        p = o -> o.r1().length() != 0 || o.r2().length() != 0;
      } else {
        p = o -> true;
      }
      try (final AsyncFastqPairWriter w = new AsyncFastqPairWriter(left, right, p)) {
        final BatchReorderingWriter<FastqPair> batchWriter = new BatchReorderingWriter<>(w);
        final Function<Batch<FastqPair>, FutureTask<?>> listRunnableFunction = batch -> new FutureTask<>(new PairAlignmentProcessor(stats, batchWriter, batch, getPairAligner()), null);
        final BatchProcessor<FastqPair> fastqPairBatchProcessor = new BatchProcessor<>(listRunnableFunction, threads, batchSize);
        fastqPairBatchProcessor.process(FastqTrim.maybeSubsample(mFlags, new FastqPairIterator(new FastqIterator(r1fq), new FastqIterator(r2fq))));
      }
      t.stop(stats.mTotal);
      t.log();
      stats.printSummary();
    }
    return 0;
  }

  private PairAligner getPairAligner() {
    final int maxReadLength = 300;
    final int seedLength = 5;
    final NgsParams ngsParams = MapParamsHelper.populateAlignerPenaltiesParams(NgsParams.builder(), mFlags).singleIndelPenalties(null).create();
    return new PairAligner(
      new UnidirectionalAdaptor(new SingleIndelSeededEditDistance(ngsParams, false, seedLength, 2, 2, maxReadLength)),
      (Integer) mFlags.getValue(MIN_OVERLAP), (Integer) mFlags.getValue(MIN_IDENTITY),
      (Integer) mFlags.getValue(LEFT_PROBE_LENGTH), (Integer) mFlags.getValue(RIGHT_PROBE_LENGTH),
      (Integer) mFlags.getValue(MIN_READ_LENGTH), mFlags.isSet(MIDPOINT_TRIM), mFlags.isSet(VERBOSE));
  }

  @Override
  public String moduleName() {
    return PETRIM_MODULE_NAME;
  }

  @Override
  public String description() {
    return "trim paired-end read FASTQ files based on read arm alignment overlap";
  }
}
