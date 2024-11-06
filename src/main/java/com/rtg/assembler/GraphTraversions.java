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

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.PathsIterator;

/**
 */
public class GraphTraversions {
  private final Map<Long, Traversion> mMap = new HashMap<>();
  GraphTraversions(Graph g) {
    for (long i = 1; i <= g.numberContigs(); ++i) {
      if (g.contigDeleted(i)) {
        continue;
      }

      final Traversion forward = new Traversion(nextContig(i, true, g), nextContig(i, false, g));
      mMap.put(i, forward);
      final HashSet<Long> reverseNext = new HashSet<>();
      for (Long l : forward.previous()) {
        reverseNext.add(-l);
      }
      final HashSet<Long> reversePrevious = new HashSet<>();
      for (Long l : forward.next()) {
        reversePrevious.add(-l);
      }
      mMap.put(-i, new Traversion(reverseNext, reversePrevious));
    }
  }

  /**
   * Retrieve a traversion for the specified contig
   * @param contigId id of the contig
   * @return the corresponding traversion
   */
  public Traversion get(long contigId) {
    return mMap.get(contigId);
  }

  static Set<Long> nextContig(long contigId, boolean forward, Graph graph) {
    final HashSet<Long> contigs = new HashSet<>();
    final PathsIterator iterator = graph.paths(contigId);
    long pathId;
    if (forward) {
      while ((pathId = iterator.nextPathId()) != 0) {
        final int index = iterator.contigIndex();
        if (graph.pathLength(pathId) > index + 1) {
          contigs.add(graph.pathContig(pathId, index + 1));
        }
      }
    } else {
      while ((pathId = iterator.nextPathId()) != 0) {
        final int index = iterator.contigIndex();
        if (index > 0) {
          contigs.add(graph.pathContig(pathId, index - 1));
        }
      }
    }
    return contigs;
  }
}
