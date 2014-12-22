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
