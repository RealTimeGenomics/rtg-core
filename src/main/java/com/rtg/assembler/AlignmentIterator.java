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

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.rtg.assembler.graph.Graph;
import com.rtg.index.hash.ExactHashFunction;

/**
 * Iterates over reads and produces alignments not an iterator because of exceptions, and laziness
 */
class AlignmentIterator {
  static class ReadAlignment {
    long mId;
    List<Set<GraphAlignment>> mFragments;

    ReadAlignment(long id, List<Set<GraphAlignment>> fragments) {
      mId = id;
      mFragments = fragments;
    }
  }
  final AsyncReadSource mReader;
  final Graph mGraph;
  final GraphAligner mAligner;
  final GraphIndex mIndex;
  final ExactHashFunction mSearchFunction;
  private final GraphMapStatistics mStatistics;
  long mId = 0;
  List<byte[]> mNextFragments;

  AlignmentIterator(AsyncReadSource reader, Graph graph, GraphAligner aligner, GraphIndex index, GraphMapStatistics statistics) {
    mReader = reader;
    mGraph = graph;
    mAligner = aligner;
    mIndex = index;
    mStatistics = statistics;
    mNextFragments = reader.nextFragments();
    mSearchFunction = mIndex.getSearchFunction();
  }


  public boolean hasNext() {
    return mNextFragments != null;
  }

  public ReadAlignment next() throws IOException {
    final List<Set<GraphAlignment>> fragmentAlignments = new ArrayList<>(mNextFragments.size());
    for (byte[] read : mNextFragments) {
      fragmentAlignments.add(getAlignments(read, mGraph, mAligner));
    }
    mNextFragments = mReader.nextFragments();
    return new ReadAlignment(mId++, fragmentAlignments);
  }

  private Set<GraphAlignment> getAlignments(byte[] read, Graph graph, GraphAligner aligner) throws IOException {
    final List<List<ContigPosition>> hits = mIndex.hits(read, graph, mSearchFunction);
    final Set<GraphAlignment> alignments = new HashSet<>();
    final int windowSize = mIndex.windowSize();
    boolean hasHits = false;
    for (int i = 0; i < read.length; ++i) {
      final List<ContigPosition> positionHits = hits.get(i);
      if (positionHits.size() > GraphMap.MAX_HITS_PER_START_POSITION) {
        continue;
      }
      nextHit: for (ContigPosition position : positionHits) {
        // Check successive windows backwards for a hit already handled
        int j = i - windowSize;
        while (j >= 0) {
          if (existsEquivalentHit(hits.get(j), i - j, position)) {
            continue nextHit;
          }
          j -= windowSize;
        }
        hasHits = true;
        final Set<GraphAlignment> aligned = aligner.align(read, i, position);
        alignments.addAll(aligned);
      }
    }
    if (!hasHits) {
      mStatistics.increment(GraphMapStatistics.Stat.NO_HITS);
    }
    return alignments;
  }

  private boolean existsEquivalentHit(final List<ContigPosition> hits, final int shift, final ContigPosition hit) {
    for (final ContigPosition h : hits) {
      if (h.mContigId == hit.mContigId && h.mPosition + shift == hit.mPosition) {
        mStatistics.increment(GraphMapStatistics.Stat.AVOIDED_ALIGNMENTS);
        return true;
      }
    }
    return false;
  }
}
