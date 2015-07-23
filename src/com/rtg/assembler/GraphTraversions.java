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
    for (long i = 1; i <= g.numberContigs(); i++) {
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
