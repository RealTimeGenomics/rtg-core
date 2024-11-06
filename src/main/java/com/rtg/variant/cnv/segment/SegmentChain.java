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
package com.rtg.variant.cnv.segment;

import java.util.ArrayList;
import java.util.TreeSet;

import com.rtg.util.diagnostic.Diagnostic;

/**
 * Chain of segments.
 */
public class SegmentChain extends ArrayList<Segment> {

  private final SegmentScorer mScorer;

  SegmentChain(final SegmentScorer scorer) {
    mScorer = scorer;
  }

  private final TreeSet<AdjacentSegments> mPriority = new TreeSet<>((a, b) -> {
    final int c = Double.compare(a.getScore(), b.getScore()); // Make first entry smallest
    if (c != 0) {
      return c;
    }
    // If the scores are equal prefer merging smaller total block length first
    final long l1 = a.getFirst().bins() + a.getSecond().bins();
    final long l2 = b.getFirst().bins() + b.getSecond().bins();
    final int c1 = Long.compare(l2, l1);
    if (c1 != 0) {
      return c1;
    }
    final int c2 = Integer.compare(a.getFirst().getStart(), b.getFirst().getStart());
    if (c2 != 0) {
      return c2;
    }
    return Integer.compare(System.identityHashCode(a), System.identityHashCode(b));
  });

  @Override
  public boolean add(final Segment segment) {
    if (!isEmpty()) {
      // Compute score to previous block
      final Segment prev = get(size() - 1);
      // Ideally prev.getEnd() <= segment.getStart(), but following still words for overlapping segments
      assert prev.getStart() < segment.getStart();
      if (segment.getSequenceName().equals(prev.getSequenceName())) {
        mPriority.add(new AdjacentSegments(mScorer.score(prev, segment), prev, segment));
      }
    }
    return super.add(segment);
  }

  void collapse() {
    if (isEmpty()) {
      return;
    }
    Diagnostic.progress("Processing: " + get(0).getSequenceName());
    double prevLambda = 0;
    double dE = Double.NEGATIVE_INFINITY;
    while (!mPriority.isEmpty()) {
      final AdjacentSegments mergeMe = mPriority.pollFirst();
      final Segment a = mergeMe.getFirst();
      final Segment b = mergeMe.getSecond();
      final int pos = indexOf(a); // this could be slow for lots of bins
      if (pos < 0 || get(pos + 1) != b) {
        continue; // Obsolete entry in priorty queue refering to already merged segments, just ignore it
      }
      // Replace segment "a" with the merged segment and delete segement "b"
      // The priority queue might still retain references to "a" and "b", but we can discard
      // those when they are encountered.
      assert get(pos + 1) == b;
      remove(b);
      final double lambda = mergeMe.getScore();
      final double delta = lambda - prevLambda;
      prevLambda = lambda;
      dE = Math.max(dE, delta);
      final Segment mergedSegment = new Segment(a, b, dE);
      set(pos, mergedSegment);
      if (pos > 0) {
        final Segment prev = get(pos - 1);
        if (mergedSegment.getSequenceName().equals(prev.getSequenceName())) {
          mPriority.add(new AdjacentSegments(mScorer.score(prev, mergedSegment), prev, mergedSegment));
        }
      }
      if (pos < size() - 1) {
        final Segment next = get(pos + 1);
        if (mergedSegment.getSequenceName().equals(next.getSequenceName())) {
          mPriority.add(new AdjacentSegments(mScorer.score(mergedSegment, next), mergedSegment, next));
        }
      }
      //System.out.println("[" + pos + "](" + a.bins() + ")(" + b.bins() + ")[" + (size() - 1) + "]  gain " + mScorer.score(a, b));
    }
  }

  @Override
  public void clear() {
    super.clear();
    mPriority.clear();
  }
}
