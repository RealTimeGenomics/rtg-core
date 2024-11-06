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

import com.rtg.assembler.graph.Contig;
import com.rtg.assembler.graph.Graph;
import com.rtg.util.IntegerOrPercentage;

/**
 */
public class GraphFlowAligner extends GraphAligner {
  static final int MISMATCH_SCORE = 1;
  static final int INDEL_SCORE = 1;

  /**
   * Construct an aligner for the provided graph where linked contigs overlap by some amount
   * @param graph the contig graph
   * @param maxMismatches the number of bases that are allowed to be different in the alignments
   * @param traverse pre calculated links between graphs
   */
  public GraphFlowAligner(Graph graph, IntegerOrPercentage maxMismatches, GraphTraversions traverse) {
    super(graph, maxMismatches, traverse);
  }

  /**
   * @param read byte representation of read
   * @param readPos position within read that the matching should start
   * @param currentContigId contig to match against
   * @param contigPosition position within the contig to start matching
   * @param increment which direction to travel along both contig and read (+1 for forwards, -1 for backwards)
   * @return a pair containing the length of the matched sequence and the number of mismatches
   */
  @Override
  ContigAlignment alignDirectedSingleContig(byte[] read, int readPos, long currentContigId, int contigPosition, int increment, int currentScore) {
    int mismatches = 0;
    final int contigLength = mGraph.contigLength(currentContigId);
    final Contig currentContig = mGraph.contig(currentContigId);
    int currentPos = readPos;
    int currentContigPosition = contigPosition;
    while (currentPos >= 0 && currentPos < read.length && currentContigPosition >= 0 && currentContigPosition < contigLength) {
//      System.err.println(forward + " contig:" + currentContigId + " -> " + currentPos + ":" + read[currentPos] + " " + currentContigPosition + ":" + currentContig.nt(currentContigPosition));
      if (read[currentPos] != currentContig.nt(currentContigPosition)) {
        boolean advanced = false;
        while (currentPos >= 0 && currentPos < read.length && read[currentPos] == read[currentPos - increment]) {
          currentPos += increment;
          mismatches += INDEL_SCORE;
          advanced = true;
        }
        while (currentContigPosition >= 0 && currentContigPosition < contigLength && currentContig.nt(currentContigPosition) == currentContig.nt(currentContigPosition - increment)) {
          currentContigPosition += increment;
          mismatches += INDEL_SCORE;
          advanced = true;
        }
        if (!advanced) {
          mismatches += MISMATCH_SCORE;
        }
        if (mismatches + currentScore > mMaxMismatches.getValue(read.length)) {
          return null;
        }
        if (advanced) {
          continue;
        }
      }
      currentPos += increment;
      currentContigPosition += increment;
    }
    return new ContigAlignment(currentPos, currentContigPosition, mismatches);
  }
}
