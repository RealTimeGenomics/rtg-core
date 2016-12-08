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

import com.rtg.alignment.ActionsHelper;
import com.rtg.alignment.GotohEditDistance;
import com.rtg.assembler.graph.MutableGraph;

/**
 */
public class PathAligner {

  private static final int GAP_OPEN_PENALTY = 1;
  private static final int GAP_EXTENSION_PENALTY = 1;
  private static final int SUBSTITUTION_PENALTY = 3;
  // todo this shift should be something like 1.2 * alignment length
  private static final int MAX_SHIFT = 500;
  static final double MAX_SHIFT_RATIO = 0.2;

  private final MutableGraph mGraph;
  private final GotohEditDistance mAligner;

  PathAligner(MutableGraph graph) {
    mGraph = graph;
    mAligner = new GotohEditDistance(GAP_OPEN_PENALTY, GAP_EXTENSION_PENALTY, SUBSTITUTION_PENALTY, 0, true);
  }

  int[] alignRight(final byte[] read, final int readStartPosition, final int readEndPosition, final byte[] template, final int templateStartPosition) {
    final int len = readEndPosition - readStartPosition;
    final int shift = maxShift(len);
    return mAligner.calculateEditDistanceFixedStart(read, readStartPosition, readEndPosition, template, templateStartPosition, Integer.MAX_VALUE, shift);
  }

  private int maxShift(int len) {
    return Math.min(MAX_SHIFT, (int) (5 + len * MAX_SHIFT_RATIO));
  }

  int[] alignLeft(final byte[] read, final int readStartPosition, final int readEndPosition, final byte[] template, final int templateEndPosition) {
    final int len = readEndPosition  - readStartPosition;
    final int shift = maxShift(len);
    return mAligner.calculateEditDistanceFixedEnd(read, readStartPosition, readEndPosition, template, templateEndPosition - len, templateEndPosition, Integer.MAX_VALUE, shift);
  }

  private int pathLength(final long pathId) {
    final int numContigs = mGraph.pathLength(pathId);
    final int overlap = mGraph.contigOverlap();
    int len = 0;
    for (int k = 0; k < numContigs; ++k) {
      len += mGraph.contigLength(mGraph.pathContig(pathId, k));
    }
    return len - (numContigs - 1) * overlap;
  }

  int align(final byte[] fragment, final int fragmentStartPosition, final long pathId, final int pathStartPosition) {
    final byte[] template = new byte[pathLength(pathId)];
    final int numContigs = mGraph.pathLength(pathId);
    final int overlap = mGraph.contigOverlap();
    for (int k = 0, j = 0, start = 0; k < numContigs; ++k) {
      final long contig = mGraph.pathContig(pathId, k);
      for (int i = start; i < mGraph.contigLength(contig); ++i) {
        template[j++] = mGraph.nt(contig, i);
      }
      start = overlap;
    }
    return ActionsHelper.alignmentScore(alignRight(fragment, fragmentStartPosition, fragment.length, template, pathStartPosition));
  }
}
