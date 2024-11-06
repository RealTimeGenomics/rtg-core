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
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

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
    GraphWriter.write(mutable.compact(), new StoreDirProxy(output), "GraphCleanup", Collections.emptySet());
  }

  static int clean(int minLength, MutableGraph mutable) {
    int totalDeleted = 0;
    int deleted = 1;
    final Set<Long> combined = new HashSet<>();
    while (deleted > 0) {
      deleted = 0;
      for (long i = 1; i <= mutable.numberContigs(); ++i) {
        if (mutable.contigDeleted(i)) {
          continue;
        }
        if (mutable.contigLength(i) > minLength) {
          continue;
        }
        final Set<Long> predecessors = MergeNodes.predecessors(mutable, i);
        final Set<Long> successors = MergeNodes.predecessors(mutable, -i);
        if (predecessors.isEmpty() || successors.isEmpty()) {
          mutable.deleteContig(i);
          ++deleted;
        }
        combined.clear();
        combined.addAll(predecessors);
        for (Long l : successors) {
          combined.add(-l);
        }
        if (combined.size() == 1) {
          for (Long l : combined) {
            if (l == i) {
              mutable.deleteContig(i);
              ++deleted;
            }
          }
        }
      }
      totalDeleted += deleted;
    }
    return totalDeleted;
  }
}
