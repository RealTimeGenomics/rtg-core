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

package com.rtg.reader;

import static com.rtg.launcher.CommonFlags.MIN_READ_LENGTH;
import static com.rtg.launcher.CommonFlags.NO_GZIP;
import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.launcher.CommonFlags.QUALITY_FLAG;
import static com.rtg.reader.FastqTrim.BATCH_SIZE;
import static com.rtg.reader.Sdf2Fasta.INTERLEAVE;
import static com.rtg.sam.SamFilterOptions.SUBSAMPLE_FLAG;
import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.cli.CommonFlagCategories.FILTERING;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.concurrent.FutureTask;
import java.util.function.Function;
import java.util.function.Predicate;

import com.rtg.alignment.SingleIndelSeededEditDistance;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.ngs.MapFlags;
import com.rtg.ngs.MapParamsHelper;
import com.rtg.ngs.NgsParams;
import com.rtg.sam.SamFilterOptions;
import com.rtg.util.Histogram;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
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
  private static final String MIDPOINT_MERGE = "midpoint-merge";

  private static final String MIN_IDENTITY = "min-identity";
  private static final String MIN_OVERLAP = "min-overlap-length";

  private static final String LEFT_PROBE_LENGTH = "left-probe-length";
  private static final String RIGHT_PROBE_LENGTH = "right-probe-length";

  private static final String DISCARD_EMPTY_READS = "discard-empty-reads";
  private static final String DISCARD_EMPTY_PAIRS = "discard-empty-pairs";

  private static final String MISMATCH_HANDLING = "mismatch-adjustment";

  private static final String MAX_READ_LENGTH = "Xmax-read-length";

  private static final String VERBOSE = "Xverbose";


  @Override
  protected void initFlags() {
    mFlags.setDescription(StringUtils.sentencify(description()));
    initFlags(mFlags);
  }

  static void initFlags(CFlags flags) {
    CommonFlagCategories.setCategories(flags);
    flags.registerRequired('l', LEFT, File.class, CommonFlags.FILE, "left input FASTQ file (AKA R1)").setCategory(INPUT_OUTPUT);
    flags.registerRequired('r', RIGHT, File.class, CommonFlags.FILE, "right input FASTQ file (AKA R2)").setCategory(INPUT_OUTPUT);
    flags.registerRequired('o', OUTPUT_FLAG, File.class, CommonFlags.FILE, "output filename prefix. Use '-' to write to standard output").setCategory(INPUT_OUTPUT);
    CommonFlags.initForce(flags);
    CommonFlags.initQualityFormatFlag(flags);
    MapFlags.initAlignerPenaltyFlags(flags);
    CommonFlags.initMinReadLength(flags);
    flags.registerOptional('L', MIN_OVERLAP, Integer.class, CommonFlags.INT, "minimum number of bases in overlap to trigger overlap trimming", 25).setCategory(SENSITIVITY_TUNING);
    flags.registerOptional('P', MIN_IDENTITY, Integer.class, CommonFlags.INT, "minimum percent identity in overlap to trigger overlap trimming", 90).setCategory(SENSITIVITY_TUNING);

    flags.registerOptional(MISMATCH_HANDLING, PairAligner.MismatchType.class, CommonFlags.STRING, "method used to alter bases/qualities at mismatches within overlap region", PairAligner.MismatchType.NONE).setCategory(FILTERING);
    flags.registerOptional('m', MIDPOINT_TRIM, "if set, trim overlapping reads to midpoint of overlap region").setCategory(FILTERING);
    flags.registerOptional('M', MIDPOINT_MERGE, "if set, merge overlapping reads at midpoint of overlap region. Result is in R1 (R2 will be empty)").setCategory(FILTERING);
    flags.registerOptional(LEFT_PROBE_LENGTH, Integer.class, CommonFlags.INT, "assume R1 starts with probes this long, and trim R2 bases that overlap into this", 0).setCategory(FILTERING);
    flags.registerOptional(RIGHT_PROBE_LENGTH, Integer.class, CommonFlags.INT, "assume R2 starts with probes this long, and trim R1 bases that overlap into this", 0).setCategory(FILTERING);
    flags.registerOptional(DISCARD_EMPTY_PAIRS, "if set, discard pairs where both reads have zero length (after any trimming)").setCategory(FILTERING);
    flags.registerOptional(DISCARD_EMPTY_READS, "if set, discard pairs where either read has zero length (after any trimming)").setCategory(FILTERING);

    flags.registerOptional(INTERLEAVE, "interleave paired data into a single output file. Default is to split to separate output files").setCategory(UTILITY);
    CommonFlags.initThreadsFlag(flags);
    SamFilterOptions.registerSubsampleFlags(flags);
    flags.registerOptional(BATCH_SIZE, Integer.class, CommonFlags.INT, "number of pairs to process per batch", 10000).setCategory(UTILITY);
    flags.registerOptional(MAX_READ_LENGTH, Integer.class, CommonFlags.INT, "maximum length of input read", 350).setCategory(FILTERING);
    flags.registerOptional(VERBOSE, "dump read alignment information to stderr").setCategory(UTILITY);
    CommonFlags.initNoGzip(flags);

    flags.setValidator(new FlagsValidator());
  }

  static class FlagsValidator implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      final File baseOutput = (File) flags.getValue(OUTPUT_FLAG);
      final BaseFile baseFile = FastqUtils.baseFile(baseOutput, !flags.isSet(NO_GZIP));
      return CommonFlags.validateInputFile(flags, LEFT)
        && CommonFlags.validateInputFile(flags, RIGHT)
        && validatedPairedOutputFiles(flags, baseOutput, baseFile)
        && flags.checkInRange(BATCH_SIZE, 1, Integer.MAX_VALUE)
        && flags.checkInRange(MIN_OVERLAP, 1, Integer.MAX_VALUE)
        && flags.checkInRange(MIN_IDENTITY, 1, 100)
        && flags.checkInRange(SUBSAMPLE_FLAG, 0.0, 1.0)
        && flags.checkInRange(MIN_READ_LENGTH, 0, Integer.MAX_VALUE)
        && flags.checkInRange(LEFT_PROBE_LENGTH, 0, Integer.MAX_VALUE)
        && flags.checkInRange(RIGHT_PROBE_LENGTH, 0, Integer.MAX_VALUE)
        && flags.checkInRange(MAX_READ_LENGTH, 1, Integer.MAX_VALUE)
        && flags.checkNand(MIDPOINT_TRIM, MIDPOINT_MERGE)
        && flags.checkNand(DISCARD_EMPTY_PAIRS, DISCARD_EMPTY_READS);
    }
  }

  static boolean validatedPairedOutputFiles(CFlags flags, File baseOutput, BaseFile baseFile) {
    if (flags.isSet(INTERLEAVE)) {
      if (!CommonFlags.validateOutputFile(flags, baseFile.file())) {
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
    return true;
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
    try (final SequenceDataSource r1fq = new FastqSequenceDataSource(FileUtils.createInputStream((File) mFlags.getValue(LEFT), true), encoding);
         final SequenceDataSource r2fq = new FastqSequenceDataSource(FileUtils.createInputStream((File) mFlags.getValue(RIGHT), true), encoding)) {

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
      final Histogram r1Lengths = new Histogram();
      final Histogram r2Lengths = new Histogram();
      try (final AsyncFastqPairWriter w = new AsyncFastqPairWriter(left, right, p, r1Lengths, r2Lengths)) {
        final BatchReorderingWriter<FastqPair> batchWriter = new BatchReorderingWriter<>(w);
        final Function<Batch<FastqPair>, FutureTask<?>> listRunnableFunction = batch -> new FutureTask<>(new PairAlignmentProcessor(stats, batchWriter, batch, getPairAligner()), null);
        final BatchProcessor<FastqPair> fastqPairBatchProcessor = new BatchProcessor<>(listRunnableFunction, threads, batchSize);
        fastqPairBatchProcessor.process(FastqTrim.maybeSubsample(mFlags, new FastqPairIterator(new FastqIterator(r1fq), new FastqIterator(r2fq))));
      }
      t.stop(stats.mTotalInput);
      t.log();
      stats.mTotalOutput = r1Lengths.sum();
      final String overlapDist = stats.mOverlapDist.getAsTsv(false, "overlap-length");
      final String fragDist = stats.mFragLengths.getAsTsv(false, "fragment-length");
      final String r1LengthDist = r1Lengths.getAsTsv(false, "read-length");
      final String r2LengthDist = r2Lengths.getAsTsv(false, "read-length");
      Diagnostic.developerLog("Overlap Distribution" + LS + overlapDist);
      Diagnostic.developerLog("Fragment Length Distribution" + LS + fragDist);
      Diagnostic.developerLog("R1 Read Length Distribution" + LS + r1LengthDist);
      Diagnostic.developerLog("R2 Read Length Distribution" + LS + r2LengthDist);
      Diagnostic.userLog(stats.printSummary());
      if (!FileUtils.isStdio(baseOutput)) {
        final BaseFile txtBase = FileUtils.getBaseFile(baseFile.getBaseFile(), false, ".txt");
        try (PrintStream summaryOut = new PrintStream(FileUtils.createTeedOutputStream(FileUtils.createOutputStream(txtBase.suffixedFile(".summary")), out))) {
          summaryOut.println(stats.printSummary());
        }
        final BaseFile tsvBase = FileUtils.getBaseFile(baseFile.getBaseFile(), false, ".tsv");
        FileUtils.stringToFile(overlapDist, tsvBase.suffixedFile(".overlap-lengths"));
        FileUtils.stringToFile(fragDist, tsvBase.suffixedFile(".fragment-lengths"));
        FileUtils.stringToFile(r1LengthDist, tsvBase.suffixedFile(".left-read-lengths"));
        FileUtils.stringToFile(r2LengthDist, tsvBase.suffixedFile(".right-read-lengths"));
      }
    }
    return 0;
  }

  private PairAligner getPairAligner() {
    final int maxReadLength = (Integer) mFlags.getValue(MAX_READ_LENGTH);
    final int seedLength = 5;
    final NgsParams ngsParams = MapParamsHelper.populateAlignerPenaltiesParams(NgsParams.builder(), mFlags).singleIndelPenalties(null).create();
    return new PairAligner(
      new SingleIndelSeededEditDistance(ngsParams, false, seedLength, 2, 2, maxReadLength),
      (Integer) mFlags.getValue(MIN_OVERLAP), (Integer) mFlags.getValue(MIN_IDENTITY),
      (Integer) mFlags.getValue(LEFT_PROBE_LENGTH), (Integer) mFlags.getValue(RIGHT_PROBE_LENGTH),
      (Integer) mFlags.getValue(MIN_READ_LENGTH),
      mFlags.isSet(MIDPOINT_TRIM),
      mFlags.isSet(MIDPOINT_MERGE),
      (PairAligner.MismatchType) mFlags.getValue(MISMATCH_HANDLING),
      mFlags.isSet(VERBOSE));
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
