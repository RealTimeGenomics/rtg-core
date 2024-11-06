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
