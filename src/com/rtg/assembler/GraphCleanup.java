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
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.UUID;

import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.io.GraphReader;
import com.rtg.assembler.graph.io.GraphWriter;
import com.rtg.util.store.StoreDirProxy;

/**
 */
public final class GraphCleanup {
  private GraphCleanup() { }
  static final int MIN_LENGTH = 200;
  static void usage() {
    System.err.println("GraphCleanup INPUT_GRAPH OUTPUT_DIRECTORY");
  }
  /**
   * Remove short isolated contigs from the input graph
   * @param args Command line flags : INPUT_GRAPH OUTPUT
   */
  public static void main(String[] args) {
    if (args.length < 2) {
      usage();
      return;
    }
    final int minLength;
    if (args.length > 2) {
      try {
        minLength = Integer.parseInt(args[2]);
      } catch (NumberFormatException e) {
        System.err.println("Invalid minimum contig length");
        usage();
        return;
      }
    } else {
      minLength = MIN_LENGTH;
    }
    final File input = new File(args[0]);
    final File output = new File(args[1]);
    if (output.exists()) {
      System.err.println("Output file already exists: " + output);
      usage();
      return;
    }
    try {
      if (!output.mkdir()) {
        System.err.println("Can't create output directory");
        usage();
        return;
      }
      process(minLength, input, output);

    } catch (IOException e) {
      System.err.println("Can't read that graph: " + e.getMessage());
      usage();
    }

  }

  private static void process(int minLength, File input, File output) throws IOException {
    final MutableGraph mutable = (MutableGraph) GraphReader.read(new StoreDirProxy(input));
    clean(minLength, mutable);
    GraphWriter.write(mutable.compact(), new StoreDirProxy(output), "GraphCleanup", Collections.<UUID>emptySet());
  }

  static int clean(int minLength, MutableGraph mutable) {
    int totalDeleted = 0;
    int deleted = 1;
    while (deleted > 0) {
      deleted = 0;
      for (long i = 1; i <= mutable.numberContigs(); i++) {
        if (mutable.contigDeleted(i)) {
          continue;
        }
        if (mutable.contigLength(i) > minLength) {
          continue;
        }
        final Set<Long> predecessors = MergeNodes.predecessors(mutable, i);
        final Set<Long> successors = MergeNodes.predecessors(mutable, -i);
        if (predecessors.size() == 0 || successors.size() == 0) {
          mutable.deleteContig(i);
          deleted++;
        }
        final Set<Long> combined = new HashSet<>();
        combined.addAll(predecessors);
        for (Long l : successors) {
          combined.add(-l);
        }
        if (combined.size() == 1) {
          for (Long l : combined) {
            if (l == i) {
              mutable.deleteContig(i);
              deleted++;
            }
          }
        }
      }
      totalDeleted += deleted;
    }
    return totalDeleted;
  }
}
