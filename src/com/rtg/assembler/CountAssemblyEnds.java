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

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.GraphFactory;
import com.rtg.assembler.graph.io.GraphReader;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.util.cli.CFlags;
import com.rtg.util.store.StoreDirProxy;

/**
 */
public class CountAssemblyEnds extends AbstractCli {

  private static final String MODULE_NAME = "countends";

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return null;
  }

  /**
   * @param args command line
   */
  public static void main(String[] args) {
    new CountAssemblyEnds().mainInit(args, System.out, System.err);
  }

  @Override
  protected void initFlags() {
    initFlagsLocal(mFlags);

  }

  protected static void initFlagsLocal(CFlags flags) {
    flags.registerRequired(File.class, CommonFlags.DIR, "input graph directory");
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream log) throws IOException {
    try (final PrintStream print = new PrintStream(out)) {
      final File in = (File) mFlags.getAnonymousValue(0);
      final Graph graph = GraphReader.read(GraphFactory.KMER, new StoreDirProxy(in));
      showEnds(print, graph);
    }
    return 0;
  }

  static void showEnds(PrintStream print, Graph graph) {
    final GraphTraversions t = new GraphTraversions(graph);
    for (long i = 1; i <= graph.numberContigs(); ++i) {
      if (graph.contigDeleted(i)) {
        continue;
      }
      final Traversion traversion = t.get(i);
      if (traversion.mNext.size() == 0) {
        print.println(i);
      }
      if (traversion.mPrevious.size() == 0) {
        print.println(-i);
      }
    }
  }
}
