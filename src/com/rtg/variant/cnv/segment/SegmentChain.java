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
