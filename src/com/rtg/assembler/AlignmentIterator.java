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
    final List<Set<GraphAlignment>> fragmentAlignments = new ArrayList<>();
    for (byte[] read : mNextFragments) {
      final Set<GraphAlignment> alignments = getAlignments(read, mGraph, mAligner);
      fragmentAlignments.add(alignments);
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
