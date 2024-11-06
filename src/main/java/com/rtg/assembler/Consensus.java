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

import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.SENSITIVITY_TUNING;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.UUID;

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.GraphFactory;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.io.GraphReader;
import com.rtg.assembler.graph.io.GraphWriter;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.LogStream;
import com.rtg.util.store.StoreDirProxy;
import com.rtg.util.store.StoreDirectory;

/**
 */
public class Consensus extends LoggedCli {

  private static final String MODULE_NAME = "consensus";
  private static final String KMER_SIZE = "kmer-size";
  /** attribute name for list of nodes that were combined to create this one */
  public static final String COMBINED = "combined";
  /** attribute description for combined */
  public static final String COMBINED_DESCRIPTION = "contigs that were combined to create this one";
  static final String CONSENSUS_READS = "consensus-reads";


  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return "reduce the complexity of a graph by removing invalid paths";
  }

  /**
   * @param args command line
   */
  public static void main(String[] args) {
    new Consensus().mainInit(args, System.out, System.err);
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(OUTPUT_FLAG);
  }

  @Override
  protected void initFlags() {
    initFlagsLocal(mFlags);
  }

  static void initCommonFlags(CFlags flags) {
    flags.registerOptional(CONSENSUS_READS, Integer.class, "Int", "number of reads necessary to form consensus along a path", 0).setCategory(SENSITIVITY_TUNING);
  }
  protected static void initFlagsLocal(CFlags flags) {
    flags.setDescription("Attempts to reduce the complexity of a graph by removing invalid paths.");
    CommonFlagCategories.setCategories(flags);
    CommonFlags.initOutputDirFlag(flags);
    initCommonFlags(flags);
    flags.registerRequired('k', KMER_SIZE, Integer.class, "Int", "size of kmer used to build the graph").setCategory(SENSITIVITY_TUNING);
    flags.registerRequired(File.class, CommonFlags.DIR, "input graph directory").setCategory(INPUT_OUTPUT);
    flags.setValidator(new Valid());
  }

  private static class Valid implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      if (!CommonFlags.validateOutputDirectory(flags)) {
        return false;
      }
      final File input = (File) flags.getAnonymousValue(0);
      if (!(input.exists() && input.isDirectory())) {
        flags.setParseMessage("Input file should be a directory in the RTG graph file format");
        return false;
      }
      return true;
    }
  }

  @Override
  protected int mainExec(OutputStream out, LogStream log) throws IOException {
    final File in = (File) mFlags.getAnonymousValue(0);
    final File output = (File) mFlags.getValue(OUTPUT_FLAG);
    final int kmerSize = (Integer) mFlags.getValue(KMER_SIZE);
    final int threshold = (Integer) mFlags.getValue(CONSENSUS_READS);
    final StoreDirProxy graphDir = new StoreDirProxy(in);
    final StoreDirProxy outputProxy = new StoreDirProxy(output);
    writeConsensus(kmerSize, threshold, graphDir, outputProxy);
    return 0;
  }

  static void writeConsensus(int kmerSize, int threshold, StoreDirectory graphDir, StoreDirectory outputProxy) throws IOException {
    final Set<UUID> sourceIds = new HashSet<>();
    sourceIds.add(GraphReader.getUUID(graphDir));
    final Map<String, String> combined = new HashMap<>();
    combined.put(COMBINED, COMBINED_DESCRIPTION);
    final Graph graph = GraphReader.read(GraphFactory.KMER, graphDir, combined, Collections.emptyMap());
    final GraphKmerAttribute mutable = (GraphKmerAttribute) graph;

    buildConsensus(kmerSize, threshold, mutable);

    GraphWriter.writeWithDeleted(mutable, outputProxy, CommandLine.getCommandLine(), sourceIds);
  }

  static void buildConsensus(int kmerSize, int threshold, GraphKmerAttribute mutable) {
    int mergeCount = 0;
    long numberContigs = 0;
    while (numberContigs != mutable.numberContigs()) {
      numberContigs = mutable.numberContigs();
      final MergeNodes merger = new MergeNodes(mutable, threshold, kmerSize - 1);
      merger.simplifyGraph();
      final ContigCollector collector = new ContigCollector(0, kmerSize, null, null, mutable);
      collector.collapse();
      ++mergeCount;
      Diagnostic.developerLog("Finished merge round: " + mergeCount);
    }
    Diagnostic.developerLog("Took " + mergeCount + " iterations to complete merging");

    if (numberContigs > 0) {
      final int[] nx = nxGraph(mutable);
      Diagnostic.info("N50: " + nx[50]);
      Diagnostic.info("N(0-100): " + Arrays.toString(nx));
    }
  }


  static int[] nxGraph(Graph g) {
    final int[] result = new int[101];
    final int totalLength = totalLength(g);
    final List<Integer> lengths = lengths(g);
    Collections.sort(lengths);
    int soFar = 0;
    for (int i = 0; i < 101; ++i) {
      result[i] = lengths.get(0);
    }
    int n = 0;
    for (int index = lengths.size() - 1 ; index >= 0 && n < 101; --index) {
      final int currentLength = lengths.get(index);
      soFar += currentLength;
      while (soFar > totalLength / 100.0 * n) {
        result[n] = currentLength;
        ++n;
      }
    }
    return result;
  }

  private static List<Integer> lengths(Graph g) {
    final List<Integer> lengths = new ArrayList<>();
    for (long i = 1; i <= g.numberContigs(); ++i) {
      if (g.contigDeleted(i)) {
        continue;
      }
      lengths.add(g.contigLength(i));
    }
    return lengths;
  }
  private static int totalLength(Graph g) {
    int totalLength = 0;
    for (Integer i : lengths(g)) {
      totalLength += i;
    }
    return totalLength;
  }
}
