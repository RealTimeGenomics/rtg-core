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

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.UUID;

import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.io.GraphWriter;
import com.rtg.launcher.NoStatistics;
import com.rtg.launcher.ParamsTask;
import com.rtg.util.Histogram;
import com.rtg.util.Pair;
import com.rtg.util.StringUtils;
import com.rtg.util.array.intindex.IntChunks;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.OneShotTimer;
import com.rtg.util.store.StoreDirProxy;

/**
 */
public class DeBruijnAssemblerTask extends ParamsTask<DeBruijnParams, NoStatistics> implements Closeable {

  /**
   * @param params       parameters for the build and search.
   * @param reportStream stream to write statistics to.
   */
  protected DeBruijnAssemblerTask(final DeBruijnParams params, final OutputStream reportStream) {
    super(params, reportStream, new NoStatistics(), null);
  }

  @Override
  public void close() {
  }
  static void writeContigs(GraphKmerAttribute graph, File outFile, Set<UUID> uuids) throws IOException {
    if (outFile.mkdir()) {
      GraphWriter.write(graph, new StoreDirProxy(outFile), CommandLine.getCommandLine(), uuids);
    } else {
      throw new IOException("Unable to create graph directory");
    }
  }


  @Override
  protected void exec() throws IOException {
    final int kmerSize = mParams.kmerSize();
    final File outputDir = mParams.directory();
    final List<ReadPairSource> sources = new ArrayList<>();
    for (final File f : mParams.inputFiles()) {
      sources.add(ReadPairSource.makeSource(f, mParams.region()));
    }

    try {
      build(kmerSize, outputDir, sources, 1);
    } finally {
      for (ReadPairSource source : sources) {
        source.close();
      }
    }
  }

  GraphKmerAttribute build(int kmerSize, File outputDir, List<ReadPairSource> sources, int numberThreads) throws IOException {
    final DeBruijnGraphBuilder dbg = new DeBruijnGraphBuilder(sources, kmerSize, mParams.useStringKmers() ? StringKmer.factory() : ByteKmer.factory(), 10, numberThreads);

    final OneShotTimer graphTimer = new OneShotTimer("Graph_build");
    final GraphKmerAttribute graph = buildGraph(dbg, kmerSize, outputDir);
    graphTimer.stopLog();

    final OneShotTimer bubbleTimer = new OneShotTimer("Graph_bubble");
    final int maxBubble = 2 * kmerSize + 10;
    Diagnostic.info("Maximum bubble length: " + maxBubble);
    Diagnostic.progress("Starting bubble popping");
    final BubbleExplorer be = new BubbleExplorer(graph, kmerSize, maxBubble, 2, mParams.mergeRatio());
    final Histogram histogram = be.popBubbles();
    Diagnostic.info("Bubble popping:" + StringUtils.LS + histogram);
    bubbleTimer.stopLog();

    final OneShotTimer outputTimer = new OneShotTimer("Graph_output");
    final HashSet<UUID> uuids = new HashSet<>();
    writeContigs(graph, new File(outputDir, "popped"), uuids);
    outputTimer.stopLog();
    final GraphKmerAttribute newGraph = new GraphKmerAttribute(kmerSize - 1, Collections.emptyMap(), Collections.emptyMap());
    graph.compact(newGraph);
    return newGraph;
  }

  private GraphKmerAttribute buildGraph(final DeBruijnGraphBuilder dbg, final int kmerSize, final File outputDir) throws IOException {
    final Set<UUID> uuids = new HashSet<>();
    final int goodThreshold = dbg.calculateGoodThreshold();
    final int minHashFrequency = mParams.minHashFrequency();
    if (minHashFrequency < 0) {
      Diagnostic.info("Using hash frequency threshold: " + goodThreshold);
      dbg.setGoodThreshold(goodThreshold);
    } else {
      Diagnostic.info("Potential hash frequency threshold: " + goodThreshold);
      Diagnostic.info("Using override hash frequency threshold: " + minHashFrequency);
      dbg.setGoodThreshold(minHashFrequency);
    }
    dbg.buildPreContigs();
    final int tipThreshold = kmerSize * 2;
    Diagnostic.info("Tip threshold: " + tipThreshold);
    writeContigs(dbg.preContigGraph(), new File(outputDir, "contigs"), uuids);
    final Pair<IntChunks, IntChunks> tipValues = dbg.calculateTipValues();

    Diagnostic.progress("Merging contigs past tips");
    final ContigCollector collector = new ContigCollector(tipThreshold, kmerSize, tipValues.getA(), tipValues.getB(), dbg.preContigGraph());
    collector.collapse();
    dbg.deleteTips(tipValues);
    Diagnostic.progress("Tips Deleted");
    writeContigs(dbg.preContigGraph(), new File(outputDir, "collapsed"), uuids);

    return dbg.preContigGraph();
  }

}
