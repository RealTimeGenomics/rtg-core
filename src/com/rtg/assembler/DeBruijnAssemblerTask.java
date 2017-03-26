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
