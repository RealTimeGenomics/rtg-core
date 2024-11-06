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
