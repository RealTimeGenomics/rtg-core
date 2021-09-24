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

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.SortedMap;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.io.GraphWriter;
import com.rtg.launcher.ParamsTask;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.store.StoreDirProxy;

/**
 */
@TestClass({"com.rtg.assembler.GraphMapTest", "com.rtg.assembler.GraphMapTaskTest"})
public class GraphMapTask extends ParamsTask<GraphMapParams, GraphMapStatistics> {
  /**
   * @param params       parameters for the build and search.
   * @param reportStream stream to write statistics to.
   */
  protected GraphMapTask(final GraphMapParams params, final OutputStream reportStream) {
    super(params, reportStream, new GraphMapStatistics(params.directory()), null);
  }

  @Override
  protected void exec() throws IOException {
    final File out = mParams.directory();
    final StoreDirProxy outProxy = new StoreDirProxy(out);
    final MutableGraph graph = mParams.graph();
    final List<ReadPairSource> reads = new ArrayList<>();
    try {
      for (File readDir : mParams.reads()) {
        final ReadPairSource readPairSource = ReadPairSource.makeSource(readDir, LongRange.NONE);
        reads.add(readPairSource);
      }
      for (File readDir : mParams.reads454()) {
        final ReadPairSource readPairSource = ReadPairSource454.makeSource(readDir, LongRange.NONE);
        reads.add(readPairSource);
      }
      for (File readDir : mParams.readsMatePair()) {
        final ReadPairSource readPairSource = ReadPairSourceMatePair.makeSource(readDir, LongRange.NONE);
        reads.add(readPairSource);
      }
      run(graph, reads, mParams, mStatistics);
      GraphWriter.write(graph, outProxy, CommandLine.getCommandLine(), Collections.emptySet());
    } finally {
      for (ReadPairSource reader : reads) {
        reader.close();
      }
    }
  }

  static void run(MutableGraph mutable, List<ReadPairSource> reads, GraphMapParams params, GraphMapStatistics stats) throws IOException {
    if (mutable.numberContigs() == 0) {
      Diagnostic.userLog("Skipping mapping because graph is empty");
      return;
    }
    final GraphIndex index = new GraphIndex(mutable, params.stepSize(), params.wordSize());
    if (params.minInsertSize() != Integer.MAX_VALUE || params.maxInsertSize() != Integer.MIN_VALUE) {
      for (ReadPairSource reader : reads) {
        reader.setMinInsertSize(params.minInsertSize());
        reader.setMaxInsertSize(params.maxInsertSize());
      }
    } else {
      GraphMap.setInsertSizes(reads, params, mutable, index);
    }

    for (ReadPairSource reader : reads) {
      reader.reset();
    }
    final ArrayList<PairConstraintWriter> writers = new ArrayList<>();
    for (int i = 0; i < params.numberThreads(); ++i) {
      writers.add(new PairConstraintWriter(new PrintStream(new FileOutputStream(params.file("constraints" + i)))));
    }
    try {
      final AsyncReadPool readPool = new AsyncReadPool("ReadForMap", reads);
      final SimpleThreadPool mapPool = new SimpleThreadPool(params.numberThreads(), "MapReads", true);
      final List<GraphMap> mapThreads = new ArrayList<>(params.numberThreads());
      for (int i = 0; i < params.numberThreads(); ++i) {
        final GraphMap mapThread = new GraphMap(index, mutable, writers.get(i), new PathTracker(new PalindromeTracker(mutable)));
        mapThreads.add(mapThread);
        mapPool.execute(new GraphMap.MapRunnable(mapThread, params.maxMismatches(), readPool.sources()));
      }
      mapPool.terminate();
      readPool.terminate();
      final List<PathTracker> trackers = new ArrayList<>(mapThreads.size());
      final List<ConstraintCache> constraints = new ArrayList<>(mapThreads.size());
      for (GraphMap thread : mapThreads) {
        stats.accumulate(thread.getStatistics());
        trackers.add(thread.getPathTracker());
        constraints.add(thread.getConstraints());
      }
      RecoverLinks.recover(constraints, mutable);
      final SortedMap<List<Long>, Integer> merged = PathTracker.merge(trackers);
      PathTracker.apply(merged, mutable);
      GraphMap.finalizeCounts(mapThreads, mutable);
    } finally {
      for (PairConstraintWriter writer : writers) {
        writer.close();
      }
    }
  }
}
