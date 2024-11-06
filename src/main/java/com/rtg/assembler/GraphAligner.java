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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.rtg.assembler.graph.Contig;
import com.rtg.assembler.graph.Graph;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.diagnostic.Diagnostic;

/**
 */
public class GraphAligner {

  private static final int SQUEEZE_MODIFIER = 0; //Integer.parseInt(System.getProperty("graphaligner.squeeze", "0"));
  static {
    Diagnostic.developerLog("Alignment squeeze modifier: " + SQUEEZE_MODIFIER);
  }
  private static final boolean FIRST_BASE_MISMATCH_BAIL = true; //Boolean.parseBoolean(System.getProperty("rtg.assembler.mismatchbail", "true"));

  final Graph mGraph;
  final GraphTraversions mTraversions;
  final IntegerOrPercentage mMaxMismatches;

  /**
   * Construct an aligner for the provided graph where linked contigs overlap by some amount
   * @param graph the contig graph
   * @param maxMismatches the number of bases that are allowed to be different in the alignments
   * @param traverse the links between contigs in the graph that you will consider
   */
  public GraphAligner(Graph graph, IntegerOrPercentage maxMismatches, GraphTraversions traverse) {
    mGraph = graph;
    mMaxMismatches = maxMismatches;
    mTraversions = traverse;
  }
  /**
   * Construct an aligner for the provided graph where linked contigs overlap by some amount
   * @param graph the contig graph
   * @param maxMismatches the number of bases that are allowed to be different in the alignments
   */
  public GraphAligner(Graph graph, IntegerOrPercentage maxMismatches) {
    this(graph, maxMismatches, new GraphTraversions(graph));
  }

  /**
   * Align the read at <code>ContigPosition</code> with the graph.
   * @param read byte representation of the read
   * @param offset offset within the read that was matched
   * @param position position within the graph that was matched
   * @return a list of the valid alignment paths
   */
  public Set<GraphAlignment> align(byte[] read, int offset, ContigPosition position) {
    final List<AlignmentChain> forwardAlignments = new ArrayList<>();
    assert read[offset] == mGraph.contig(position.mContigId).nt(position.mPosition);
    alignDirectedContigs(read, offset, position.mContigId, position.mPosition, true, null, forwardAlignments);
    if (forwardAlignments.isEmpty()) {
      return Collections.emptySet();
    }
    final List<AlignmentChain> reverseAlignments = new ArrayList<>();
    alignDirectedContigs(read, offset, position.mContigId, position.mPosition, false, null, reverseAlignments);
    if (reverseAlignments.isEmpty()) {
      return Collections.emptySet();
    }
    return squeezePaths(joinPaths(forwardAlignments, reverseAlignments, read.length), mGraph);
  }
  Set<GraphAlignment> squeezePaths(Set<GraphAlignment> initial, Graph g) {
    final Set<GraphAlignment> result = new HashSet<>(initial.size());
    for (GraphAlignment alignment : initial) {
      result.add(sqeezePath(alignment, g));
    }
    return result;
  }

  GraphAlignment sqeezePath(GraphAlignment alignment, Graph g) {
    GraphAlignment result = alignment;
    int startLength = g.contigLength(result.contigs().get(0));
    final int contigOverlap = getContigOverlap();
    while (result.contigs().size() > 1 && result.startPosition() >= startLength - contigOverlap - SQUEEZE_MODIFIER) {
      final int newStart = result.startPosition() - startLength + contigOverlap;
      final List<Long> newContigs = result.contigs().subList(1, result.contigs().size());
      result = new GraphAlignment(newStart, result.endPosition(), newContigs, result.mScore, g);
      startLength = g.contigLength(result.contigs().get(0));
    }
    while (result.contigs().size() > 1 && result.endPosition() < contigOverlap + SQUEEZE_MODIFIER) {
      final int previousLength = g.contigLength(result.contigs().get(result.contigs().size() - 2));
      final int newEnd = result.endPosition() + previousLength - contigOverlap;
      final List<Long> newContigs = result.contigs().subList(0, result.contigs().size() - 1);
      result = new GraphAlignment(result.startPosition(), newEnd, newContigs, result.mScore, g);
    }


    return result;
  }

  private Set<GraphAlignment> joinPaths(List<AlignmentChain> forward, List<AlignmentChain> reverse, int readLength) {
    final Set<GraphAlignment> paths = new HashSet<>();
    for (AlignmentChain aForward : forward) {
      for (AlignmentChain aReverse : reverse) {
        if (aForward.mScore + aReverse.mScore <= mMaxMismatches.getValue(readLength)) {
          paths.add(new GraphAlignment(aForward, aReverse, mGraph));
        }
      }
    }
    return paths;
  }

  private void alignDirectedContigs(byte[] read, int startReadPosition, long startContig, int startContigPosition, boolean forward, AlignmentChain soFar, List<AlignmentChain> results) {
    final int increment = forward ? 1 : -1;
    final ContigAlignment contigAlignment = alignDirectedSingleContig(read, startReadPosition, startContig, startContigPosition, increment, soFar == null ? 0 : soFar.mScore);
    if (contigAlignment == null) {
      return;
    }
//    final int alignLength = Math.abs(startReadPosition - contigAlignment.mReadEnd);
    final int mismatches = contigAlignment.mMismatches;

    final int contigAlignEnd = contigAlignment.mContigEnd - increment;
    final int readAlignEnd = contigAlignment.mReadEnd - increment; //startReadPosition + increment * (alignLength - 1);
    final AlignmentSection section;
    if (forward) {
      section = new AlignmentSection(startContig, startContigPosition, contigAlignEnd);
    } else {
      section = new AlignmentSection(startContig, contigAlignEnd, startContigPosition);
    }

    final AlignmentChain nextLink = new AlignmentChain((soFar == null ? 0 : soFar.mScore) + mismatches, section, soFar);
    if (readAlignEnd == 0 || readAlignEnd == read.length - 1) {
      results.add(nextLink);
      return;
    }
    final Traversion next = mTraversions.get(startContig);
    for (long contigId : forward ? next.next() : next.previous()) {
      final int nextAlignStart = forward ? getContigOverlap() : mGraph.contigLength(contigId) - getContigOverlap() - 1;
      alignDirectedContigs(read, readAlignEnd + increment, contigId, nextAlignStart, forward, nextLink, results);
    }
  }

  public int getContigOverlap() {
    return mGraph.contigOverlap();
  }

  static class ContigAlignment {
    int mMismatches;
    int mReadEnd;
    int mContigEnd;
    ContigAlignment(int readEnd, int contigEnd, int mismatches) {
      mMismatches = mismatches;
      mReadEnd = readEnd;
      mContigEnd = contigEnd;
    }

    @Override
    public String toString() {
      return "ContigAlignment{" + "mMismatches=" + mMismatches + ", mReadEnd=" + mReadEnd + ", mContigEnd=" + mContigEnd + '}';
    }
  }

  /**
   * @param read byte representation of read
   * @param readPos position within read that the matching should start
   * @param currentContigId contig to match against
   * @param contigPosition position within the contig to start matching
   * @param increment which direction to travel along both contig and read (+1 for forwards, -1 for backwards)
   * @param currentScore contain the score computed so far
   * @return a pair containing the length of the matched sequence and the number of mismatches
   */
  ContigAlignment alignDirectedSingleContig(byte[] read, int readPos, long currentContigId, int contigPosition, int increment, int currentScore) {
    int mismatches = 0;
    final int contigLength = mGraph.contigLength(currentContigId);
    final Contig currentContig = mGraph.contig(currentContigId);
    int currentPos = readPos;
    int currentContigPosition = contigPosition;
    while (currentPos >= 0 && currentPos < read.length && currentContigPosition >= 0 && currentContigPosition < contigLength) {

      if (read[currentPos] != currentContig.nt(currentContigPosition)) {
        if (FIRST_BASE_MISMATCH_BAIL && contigPosition == currentContigPosition) {
          return null;
        }
        ++mismatches;
        if (mismatches + currentScore > mMaxMismatches.getValue(read.length)) {
          return null;
        }
      }
      currentPos += increment;
      currentContigPosition += increment;
    }
    return new ContigAlignment(currentPos, currentContigPosition, mismatches);
  }
}
