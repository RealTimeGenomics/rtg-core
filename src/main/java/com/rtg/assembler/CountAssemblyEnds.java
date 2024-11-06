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
      if (traversion.mNext.isEmpty()) {
        print.println(i);
      }
      if (traversion.mPrevious.isEmpty()) {
        print.println(-i);
      }
    }
  }
}
