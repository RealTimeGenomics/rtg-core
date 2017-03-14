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
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.List;

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.assembler.graph.io.GraphReader;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.util.cli.CFlags;
import com.rtg.util.store.StoreDirProxy;

/**
 */
public class ExtractPath extends AbstractCli {

  private static final String PATH = "path";
  private static final String KMER_SIZE = "kmer-size";
  private static final String START = "start";
  private static final String END = "end";


  @Override
  public String moduleName() {
    return "ExtractPath";
  }

  @Override
  public String description() {
    return null;
  }

  @Override
  protected void initFlags() {
    initFlagsLocal(mFlags);
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final File graphDir = (File) mFlags.getAnonymousValue(0);
    final Graph graph = GraphReader.read(new StoreDirProxy(graphDir));
    final long pathId = (Long) mFlags.getValue(PATH);
    final int kmer = (Integer) mFlags.getValue(KMER_SIZE);
    if (graph.numberPaths() <= pathId) {
      throw new IllegalArgumentException("that path doesn't exist");
    }
    final int start;
    if (mFlags.isSet(START)) {
      start = (Integer) mFlags.getValue(START);
    } else {
      start = -1;
    }
    final int end;
    if (mFlags.isSet(END)) {
      end = (Integer) mFlags.getValue(END);
    } else {
      end = -1;
    }
    run(out, graph, pathId, kmer, start, end);
    return 0;
  }

  static void run(OutputStream out, Graph graph, long pathId, int kmer, int start, int end) {
    final StringBuilder sb = new StringBuilder();
    final List<Long> path = MergeNodes.getPath(graph, pathId);
    sb.append(ContigString.contigSequenceString(graph.contig(path.get(0))).substring(0, kmer - 1));
    for (long contigId : path) {
      sb.append(ContigString.contigSequenceString(graph.contig(contigId)).substring(kmer - 1));
      sb.append("|");
    }
    // remove trailing "|"
    sb.setLength(sb.length() - 1);
    if (start > -1) {
      if (start < graph.contigLength(path.get(0))) {
        sb.delete(0, start);
      }
    }
    if (end > -1) {
      final int contigLength = graph.contigLength(path.get(path.size() - 1));
      if (end < contigLength) {
        final int tail = contigLength - end - 1;
        sb.setLength(sb.length() - tail);
      }
    }
    try (final PrintStream output = new PrintStream(out)) {
      output.println(sb);
    }
  }

  protected static void initFlagsLocal(CFlags flags) {
    flags.registerExtendedHelp();
    flags.setDescription("Outputs the sequence pertaining to the specified path.");
    flags.registerRequired(File.class, CommonFlags.DIR, "input graph directory");
    flags.registerRequired('p', PATH, Long.class, CommonFlags.INT, "path that will be traced through the graph and output");
    flags.registerRequired('k', KMER_SIZE, Integer.class, CommonFlags.INT, "size of the kmer the graph was constructed with");
    flags.registerOptional('s', START, Integer.class, CommonFlags.INT, "trim the first sequence up to this position");
    flags.registerOptional('e', END, Integer.class, CommonFlags.INT, "trim the final sequence at this position");
  }

  /**
   * @param args command line
   */
  public static void main(String[] args) {
    new ExtractPath().mainInit(args, System.out, System.err);
  }
}
