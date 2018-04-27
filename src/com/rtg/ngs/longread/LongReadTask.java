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
package com.rtg.ngs.longread;

import java.io.IOException;

import com.rtg.index.Add;
import com.rtg.index.Finder;
import com.rtg.index.Index;
import com.rtg.index.IndexFilterMethod;
import com.rtg.index.IndexUtils;
import com.rtg.index.hash.ExactHashFunction;
import com.rtg.index.hash.HashFunction;
import com.rtg.index.hash.HashLoop;
import com.rtg.index.hash.InExactHashFunction;
import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.index.params.CreateParams;
import com.rtg.index.queue.IndexQueues;
import com.rtg.launcher.BuildParams;
import com.rtg.launcher.HashingRegion;
import com.rtg.position.FinderPositionOutput;
import com.rtg.position.PositionUtils;
import com.rtg.position.PositionWriterFactory;
import com.rtg.position.SearchIncrementalHashLoop;
import com.rtg.position.SearchResetHashLoop;
import com.rtg.position.output.GapBucketsInfo;
import com.rtg.position.output.OutputFormatType;
import com.rtg.position.output.OutputProcessorWrapper;
import com.rtg.position.output.PositionOutput;
import com.rtg.position.output.PositionOutputParams;
import com.rtg.position.output.PositionParams;
import com.rtg.position.output.PositionWriter;
import com.rtg.reader.SequencesReader;
import com.rtg.reference.Sex;
import com.rtg.usage.UsageMetric;
import com.rtg.util.IORunnable;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.SizeSplit;
import com.rtg.util.StringUtils;
import com.rtg.util.array.ImmutableIntArray;
import com.rtg.util.array.SingleValueIntArray;
import com.rtg.util.array.WrappedIntArray;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.OneShotTimer;

/**
 * Takes subjects, queries and writes to "out" in output directory.
 * Also accepts the various other parameters required.
 * Specialized version for long reads only called from <code>RamMap via NgsTask</code>.
 */
public final class LongReadTask {

  //  private static final boolean SIMPLE_TUBE = false; //Boolean.valueOf(System.getProperty("simple.tube", "false"));

  private LongReadTask() { }

  /**
   * @param params parameters for the build and search.
   * @param output where grouped hits are sent for further processing
   * @param usageMetric accumulates the count of nucleotides read.
   * @throws IOException if an I/O error occurs
   */
  static void execute(final PositionParams params, final OutputProcessor output, final UsageMetric usageMetric) throws IOException {
    final OneShotTimer fullTimer = new OneShotTimer("total_time");
    final Index index = build(params, usageMetric, params.ngsParams().indexFilter());
    search(params, output, index);
    fullTimer.stopLog();
    Diagnostic.userLog("Index search performance " + StringUtils.LS + index.perfString());
  }

  /**
   * Search an index with given parameters
   * @param params search parameters
   * @param output output processor
   * @param index index to search
   * @throws IOException if an I/O error occurs
   */
  public static void search(final PositionParams params, final OutputProcessor output, final Index index) throws IOException {
    //Search the queries and write hits
    //open the query reader
    final OneShotTimer writeTimer = new OneShotTimer("LR_BS_search");
    final ImmutableIntArray readLengths;
    final SequencesReader reader = params.build().sequences().reader();


    if (params.ngsParams().paired()) {
      final SequencesReader readerSecond = params.buildSecond().sequences().reader();
      if (reader.minLength() == reader.maxLength() && reader.minLength() == readerSecond.minLength() && readerSecond.minLength() == readerSecond.maxLength()) {
        assert (reader.numberSequences() + readerSecond.numberSequences()) <= Integer.MAX_VALUE;
        readLengths = new SingleValueIntArray((int) reader.minLength(), (int) (reader.numberSequences() + readerSecond.numberSequences()));
      } else {
        final int[] seqLengths = reader.sequenceLengths(0, reader.numberSequences());
        final int[] seqLengthsSecond = readerSecond.sequenceLengths(0, readerSecond.numberSequences());
        readLengths = new WrappedIntArray(seqLengths, seqLengthsSecond);
      }
    } else {
      final int[] seqLengths = reader.sequenceLengths(0, reader.numberSequences());
      readLengths = new WrappedIntArray(seqLengths);
    }
    search(params, index, readLengths, output);
    writeTimer.stopLog();
  }

  /**
   * Build index using given parameters
   * @param params index parameters
   * @param usageMetric accumulates the count of nucleotides read.
   * @param hashFilter handle used to filter hashes from index
   * @return index
   * @throws IOException if an I/O error occurs
   */
  public static Index build(final PositionParams params, final UsageMetric usageMetric, IndexFilterMethod hashFilter) throws IOException {
    Diagnostic.developerLog("Estimated usage of memory" + StringUtils.LS + PositionUtils.memToString(params));
    final OneShotTimer buildTimer = new OneShotTimer("LR_BS_build");
    final int numberThreads = params.numberThreads();
    final OneShotTimer queueTimer = new OneShotTimer("LR_BS_build_queue");
    final SimpleThreadPool pool = new SimpleThreadPool(numberThreads, "Build", true);
    pool.enableBasicProgress(numberThreads);
    final SizeSplit ss = new SizeSplit((int) params.build().sequences().reader().numberSequences(), numberThreads);
    final CreateParams indexParams;
    final boolean paired = params.ngsParams().paired();
    indexParams = params.indexParams();
    final int ipBits = indexParams.initialPointerBits();
    final IndexQueues queues = new IndexQueues(numberThreads, indexParams.hashBits(), indexParams.size(), indexParams.valueBits(), ipBits);
    for (int i = 0; i < numberThreads; ++i) {
      final int start = ss.start(i);
      final int end = ss.start(i + 1);
      final HashingRegion region = new HashingRegion(start, end);

      final BTask buildTask;
      if (paired) {
        buildTask = new BuildPairedTask(params.subBuild(region), queues.queue(i), usageMetric);
      } else {
        buildTask = new BuildTask(params.subBuild(region), queues.queue(i), usageMetric);
      }
      pool.execute(buildTask);
    }
    pool.terminate();
    queueTimer.stopLog();
    final OneShotTimer freezeTimer = new OneShotTimer("LR_BS_freeze");
    final Index index = IndexUtils.createIndex(indexParams, hashFilter, numberThreads);
    queues.freeze(index);
    freezeTimer.stopLog();
    buildTimer.stopLog();
    //System.err.println(index);
    Diagnostic.userLog("Index statistics " + StringUtils.LS + index.infoString());

    return index;
  }

  private abstract static class BTask implements IORunnable {
    protected final PositionParams mParams;
    protected final Add mQueue;
    protected final UsageMetric mUsageMetric;

    BTask(final PositionParams params, final Add queue, final UsageMetric usageMetric) {
      mParams = params;
      mQueue = queue;
      mUsageMetric = usageMetric;
    }
  }

  private static class BuildTask extends BTask {

    /**
     * @param params parameters.
     * @param queue where hash/id pairs are added.
     * @param usageMetric counter for accumulating number of nucleotides read.
     */
    BuildTask(final PositionParams params, final Add queue, final UsageMetric usageMetric) {
      super(params, queue, usageMetric);
    }

    @Override
    public void run() throws IOException {
      final OneShotTimer timer = new OneShotTimer("LR_BS_build_queue_" + Thread.currentThread().getName().replace(" ", "_"));
      final HashLoop hl = PositionUtils.makeBuild(mQueue, mParams.build(), mParams.indexParams().windowBits());
      final byte[] buffer = PositionUtils.makeBuffer(mParams.build().sequences().maxLength());
      mUsageMetric.incrementMetric(hl.execLoop(mParams.build().sequences(), buffer));
      timer.stopLog();
    }

  }

  private static class BuildPairedTask extends BTask {

    /**
     * @param params parameters.
     * @param queue where hash/id pairs are added.
     * @param usageMetric counter for accumulating number of nucleotides read.
     */
    BuildPairedTask(final PositionParams params, final Add queue, final UsageMetric usageMetric) {
      super(params, queue, usageMetric);
    }

    @Override
    public void run() throws IOException {
      final OneShotTimer timer = new OneShotTimer("LR_BS_build_queue_" + Thread.currentThread().getName().replace(" ", "_"));
      final BuildParams buildFirst = mParams.build();
      final BuildParams buildSecond = mParams.buildSecond();
      final HashLoop hl = PositionUtils.makePairedBuild(mQueue, buildFirst, buildSecond, mParams.indexParams().windowBits());
      final long lLeft = mParams.build().sequences().maxLength();
      final long lRight = mParams.buildSecond().sequences().maxLength();
      final byte[] buffer = PositionUtils.makeBuffer(Math.max(lLeft, lRight));

      hl.setSide(false);
      final long length1 = hl.execLoop(mParams.build().sequences(), buffer);
      hl.setSide(true);
      final long length2 = hl.execLoop(mParams.buildSecond().sequences(), buffer);
      mUsageMetric.incrementMetric(length1 + length2);
      timer.stopLog();
    }
  }

  private static HashingRegion[] ngsRegions(long length, int nt, int multiplier, HashingRegion superRegion, SequencesReader sr, long padding, Sex sex) throws IOException {
    final long rStart = superRegion.getStart() == HashingRegion.MISSING ? 0 : superRegion.getStart();
    final long rEnd = rStart + length;
    return HashingRegion.splitWorkload(sr, sex, rStart, rEnd, nt * multiplier, HashingRegion.DEFAULT_MIN_CHUNK_SIZE, padding);
  }

  private static void search(PositionParams params, Index index, ImmutableIntArray readLengths, OutputProcessor output) throws IOException {
    final HashingRegion region = params.search().sequences().region();
    final long length;
    if (region != HashingRegion.NONE) {
      length = region.getEnd() - region.getStart();
    } else {
      length = params.search().sequences().reader().numberSequences();
    }
    final int nt = params.numberThreads();
    assert nt >= 1;
    final SimpleThreadPool pool = new SimpleThreadPool(nt, "LongReadSearch", true);
    final int multiplier;
    if (nt == 1) {
      multiplier = 1;
    } else {
      multiplier = params.ngsParams().threadMultiplier();
    }
    final Sex sex = params.ngsParams().sex();
    final HashingRegion[] regions = ngsRegions(length, nt, multiplier, region, params.search().sequences().reader(), params.ngsParams().calculateThreadPadding(), sex);

    //System.err.println(params + "=params " + params.output() + "=params.output() " + params.output().maxGap() + "=params.output().maxGap() ");
    //final GapBucketsInfo bucketInfoReverse = params.output().maxGap() == null ? null : new GapBucketsInfo(params, 1, false, readLengths); //it seems it wasn't necessary
    //System.err.println(" start=" + start + " end=" + end + " length=" + length + " n=" + n + " m=" + m);
    final PositionWriterFactory writerFactory = new PositionWriterFactory.NgsPositionWriterFactory(output);
    final GapBucketsInfo bucketInfo;
    if (params.ngsParams().paired()) {
      bucketInfo = params.output().maxGap() == null ? null : new GapBucketsInfo(params, 1, true);
    } else {
      bucketInfo = params.output().maxGap() == null ? null : new GapBucketsInfo(params, 1);
    }
    int count = 0;
    pool.enableBasicProgress(regions.length);
    for (final HashingRegion r : regions) {
      final PositionParams subParams = params.subSearch(r);
      Diagnostic.developerLog("LongReadTask Thread create: " + r);
      pool.execute(new SubSearch(index, readLengths, count++, subParams, writerFactory, bucketInfo));
    }
    pool.terminate();
  }

  private static class SubSearch implements IORunnable {

    private static class PositionFinder extends Finder {

      /** The maximum number of hit positions in a build sequence (based on max sequence length and step size) */
      private final int mMxs;
      private final PositionOutput mOutput;
      private final int mStepSize;
      private final boolean mReverseMode;
      private final ImmutableIntArray mBuildLengths;
      private final int mWindowLength;

      PositionFinder(final int mxs, final PositionOutput output, final int stepSize, final boolean reverseMode, final ImmutableIntArray buildLengths, final int windowLength) {
        mMxs = mxs;
        mOutput = output;
        mStepSize = stepSize;
        mReverseMode = reverseMode;
        mBuildLengths = buildLengths;
        mWindowLength = windowLength;
      }

      @Override
      public boolean found(final long id) throws IOException {
        //System.err.println("found id=" + id);
        assert mMxs > 0; //if we actually get here then this will be true
        final int seqId = (int) (id / mMxs);
        final int stepPosition = (int) (id % mMxs);
        assert stepPosition >= 0 && stepPosition < mMxs; // : "id=" + id + " mMxs=" + mMxs + " mStepSize=" + mStepSize;

        final int position = stepPosition * mStepSize;
        final int posn;
        if (mReverseMode) {
          posn = mBuildLengths.get(seqId) - position - mWindowLength;
        } else {
          posn = position;
        }
        mOutput.hit(seqId, posn);
        return true;
      }
    }

    private final Index mIndex;
    private final ImmutableIntArray mReadLengths;
    private final int mT;
    private final PositionParams mSubParams;
    private final PositionWriterFactory mWriterFactory;
    private final GapBucketsInfo mBucketInfo;

    SubSearch(final Index index, final ImmutableIntArray readLengths, final int t, final PositionParams subParams, final PositionWriterFactory writerFactory, final GapBucketsInfo bucketInfo) {
      mIndex = index;
      mReadLengths = readLengths;
      mT = t;
      mSubParams = subParams;
      mWriterFactory = writerFactory;
      mBucketInfo = bucketInfo;
    }

    @Override
    public void run() throws IOException {
      Diagnostic.developerLog("LongReadTask Thread start: " + mSubParams.search().sequences().region());
      final OneShotTimer timer = new OneShotTimer("Search_" + mT);
      try {
        final BuildParams buildParams = mSubParams.build();
        final int mxs;
        if (mSubParams.ngsParams().paired()) {
          mxs = Math.max(PositionUtils.maxMatches(buildParams), PositionUtils.maxMatches(mSubParams.buildSecond()));
        } else {
          mxs = PositionUtils.maxMatches(buildParams);
        }
        assert mxs >= 0; //if zero then no actual hits.
        final int stepSize = buildParams.stepSize();
        final PositionWriter writer;
        final PositionOutput output;
        if (mSubParams.ngsParams().paired()) {
          writer = mWriterFactory.makeNgs(mSubParams.ngsParams().buildFirstParams(), mSubParams.ngsParams().buildSecondParams(),
              mSubParams.search().sequences(), mSubParams.search().sequences().region());
        } else {
          writer = mWriterFactory.makeNgs(mSubParams.build().sequences(), null,
              mSubParams.search().sequences(), mSubParams.search().sequences().region());
        }
        final PositionOutputParams op = mSubParams.output();
        final OutputFormatType f = op.format();
        output = f.output(mSubParams, mBucketInfo, mReadLengths, System.out, null, writer, mIndex.maxHashCount());

        final FinderPositionOutput outputVars = new FinderPositionOutput(new PositionFinder(mxs, output, stepSize, false, null, buildParams.windowSize()), output);

        final PositionOutput temp = output.reverseClone();
        final FinderPositionOutput outputVarsReverse = new FinderPositionOutput(new PositionFinder(mxs, temp, stepSize, true, mReadLengths, buildParams.windowSize()), temp);

        final int winBits = mSubParams.indexParams().windowBits();
        final HashLoop queryHashLoop;
        if (winBits > 64) {
          final HashFunction function = new InExactHashFunction(buildParams.windowSize());
          queryHashLoop = new SearchResetHashLoop(mSubParams.search().windowSize(), mSubParams.search().stepSize(), function, outputVars, outputVarsReverse, mIndex, true);
        } else {
          final ExactHashFunction function = new ExactHashFunction(buildParams, true);
          // xxx SpecializedIncrementalHashLoop claims to be optimized, but does it really help?
          //queryHashLoop = new SpecializedIncrementalHashLoop(mSubParams.search().stepSize(), function, outputVars, outputVarsReverse,  mIndex, true);
          queryHashLoop = new SearchIncrementalHashLoop(mSubParams.search().stepSize(), function, outputVars, outputVarsReverse,  mIndex, true);
        }
        final int padding = mSubParams.ngsParams().calculateThreadPadding(); // Amount added either side
        final byte[] buffer = new byte[(int) mSubParams.search().sequences().region().longestSubSequence(mSubParams.search().sequences().reader()) + 2 * padding];
        //System.err.println(mSubParams + "=mSubParams");
        queryHashLoop.setThreadPadding(padding);
        queryHashLoop.execLoop(mSubParams.search().sequences(), buffer);
        ((OutputProcessorWrapper) writer).getEnclosed().threadFinish();
        Diagnostic.developerLog("LongReadTask Thread finish: " + mSubParams.search().sequences().region());
      } finally {
        mSubParams.close();
      }
      timer.stopLog();
    }
  }
}

