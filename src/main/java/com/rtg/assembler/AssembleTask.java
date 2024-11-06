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
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.GraphFactory;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.io.GraphReader;
import com.rtg.assembler.graph.io.GraphWriter;
import com.rtg.launcher.ParamsTask;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.store.StoreDirProxy;

/**
 */
final class AssembleTask extends ParamsTask<AssembleParams, GraphMapStatistics> {
  /**
   * @param params       parameters for the build and search.
   * @param reportStream stream to write statistics to.
   * @param stats        statistics object to populate.
   */
  protected AssembleTask(final AssembleParams params, final OutputStream reportStream, GraphMapStatistics stats) {
    super(params, reportStream, stats, null);
  }

  private static LongRange adjustRegionRestriction(final LongRange region, final File f) throws IOException {
    if (ReaderUtils.isPairedEndDirectory(f)) {
      return SequencesReaderFactory.resolveRange(ReaderUtils.getLeftEnd(f), region);
    } else {
      return SequencesReaderFactory.resolveRange(f, region);
    }
  }

  /**
   * Run all steps in assemble process
   * @param assembleParams parameters containing all the assembly options
   * @param printStream output goes here
   * @param stats statistics object to populate.
   * @throws IOException when some terrible tragedy occurs in the assembly IO
   */
  public static void assemble(AssembleParams assembleParams, PrintStream printStream, GraphMapStatistics stats) throws IOException {
    // TODO region handling is kind of weird here, given potential for more than one read set
    // currently this assumes any restriction applies to all read sets.
    final LongRange region = assembleParams.region();
    final List<File> reads = assembleParams.reads();
    final List<File> reads454 = assembleParams.reads454();
    final List<File> readsMatePair = assembleParams.readsMatePair();
    final List<ReadPairSource> readPairSources = new ArrayList<>(reads.size() + reads454.size() + readsMatePair.size());
    for (File f : reads) {
      final LongRange r = adjustRegionRestriction(region, f);
      readPairSources.add(ReadPairSource.makeSource(f, r));
    }
    for (File f : reads454) {
      final LongRange r = adjustRegionRestriction(region, f);
      readPairSources.add(ReadPairSource454.makeSource(f, r));
    }
    for (File f : readsMatePair) {
      final LongRange r = adjustRegionRestriction(region, f);
      readPairSources.add(ReadPairSourceMatePair.makeSource(f, r));
    }
    try {
      final GraphKmerAttribute graph = makeGraph(assembleParams, printStream, readPairSources);
      mapReads(assembleParams, stats, readPairSources, graph);
      writeGraph(assembleParams.directory(), graph, "mapped", false);
      Diagnostic.userLog("Starting improve single");
      if (assembleParams.minReadCount() > -1) {
        FilterPaths.improveSingle(graph, assembleParams.minReadCount());
      }
      Diagnostic.userLog("Starting improve multiple");
      if (assembleParams.minPathReads() > -1) {
        FilterPaths.improveMultiple(graph, assembleParams.minPathReads());
      }
      Diagnostic.userLog("Starting consensus");
      Consensus.buildConsensus(assembleParams.kmerSize(), assembleParams.consensusThreshold(), graph);

      final Graph sorted = GraphSorter.sortedGraph(graph);
      writeGraph(assembleParams.directory(), sorted, "final", false);

      writeGraph(assembleParams.directory(), graph, "unfiltered_final", true);
    } finally {
      for (ReadPairSource source : readPairSources) {
        source.close();
      }
    }
  }
  static void writeGraph(File directory, Graph g, String fileName, boolean withDeleted) throws IOException {
    final File deletes = new File(directory, fileName);
    if (!deletes.mkdir()) {
      throw new NoTalkbackSlimException(ErrorType.DIRECTORY_NOT_CREATED, deletes.getPath());
    }
    if (withDeleted) {
      GraphWriter.writeWithDeleted(g, new StoreDirProxy(deletes), CommandLine.getCommandLine(), Collections.emptySet());
    } else {
      GraphWriter.write(g, new StoreDirProxy(deletes), CommandLine.getCommandLine(), Collections.emptySet());
    }

  }

  private static void mapReads(AssembleParams assembleParams, GraphMapStatistics stats, List<ReadPairSource> readPairSources, GraphKmerAttribute graph) throws IOException {
    if (readPairSources.size() > 0) {
      final GraphMapParams params = GraphMapParams.builder()
          .wordSize(assembleParams.wordSize())
          .stepSize(assembleParams.stepSize())
          .maxMismatches(assembleParams.maxMismatches())
          .minInsertSize(assembleParams.minInsertSize())
          .maxInsertSize(assembleParams.maxInsertSize())
          .numberThreads(assembleParams.numberThreads())
          .directory(assembleParams.directory())
          .create();
      GraphMapTask.run(graph, readPairSources, params, stats);
    }
  }

  private static GraphKmerAttribute makeGraph(AssembleParams assembleParams, PrintStream printStream, List<ReadPairSource> readPairSources) throws IOException {
    final GraphKmerAttribute graph;
    if (assembleParams.graph() == null) {
      final DeBruijnParams params = DeBruijnParams.builder().mergeRatio(assembleParams.mergeRatio()).minHashFrequency(assembleParams.minHashFrequency()).region(assembleParams.region()).create();
      final File buildDirectory = assembleParams.file("build");
      if (!buildDirectory.mkdir()) {
        throw new NoTalkbackSlimException(ErrorType.DIRECTORY_NOT_CREATED, buildDirectory.getPath());
      }
      graph = new DeBruijnAssemblerTask(params, printStream).build(assembleParams.kmerSize(), buildDirectory, readPairSources, assembleParams.numberThreads());
    } else {
      graph = (GraphKmerAttribute) GraphReader.read(GraphFactory.KMER, new StoreDirProxy(assembleParams.graph()));
    }
    graph.addContigAttribute(GraphKmerAttribute.READ_COUNT, GraphKmerAttribute.READ_COUNT_DESCRIPTION);
    graph.addContigAttribute(Consensus.COMBINED, Consensus.COMBINED_DESCRIPTION);
    graph.addPathAttribute(GraphKmerAttribute.READ_COUNT, GraphKmerAttribute.READ_COUNT_DESCRIPTION);
    return graph;
  }

  @Override
  protected void exec() throws IOException {
    try (final PrintStream print = new PrintStream(mReportStream)) {
      assemble(mParams, print, getStatistics());
    }
  }
}
