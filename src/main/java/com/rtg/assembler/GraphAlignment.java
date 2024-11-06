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
import java.util.List;

import com.rtg.assembler.graph.Contig;
import com.rtg.assembler.graph.Graph;
import com.rtg.util.LongUtils;
import com.rtg.util.Utils;

/**
 */
public class GraphAlignment {
  final List<Long> mContigs;
  final int mScore;
  final int mStartPosition;
  final int mEndPosition;
  final Graph mGraph;

  @Override
  public String toString() {
    return "GraphAlignment{" + "mContigs=" + mContigs + ", mStartPosition=" + mStartPosition + ", mEndPosition=" + mEndPosition +  ", mScore=" + mScore + '}';
  }

  /**
   * Represents the path through a graph that an alignment travels.
   * @param parts the <code>AlignmentSection</code>s containing the contigs that have been aligned
   * @param score the alignment score
   * @param graph the graph this alignment is against
   */
  public GraphAlignment(List<AlignmentSection> parts, int score, Graph graph) {
    mScore = score;
    mGraph = graph;
    mStartPosition = parts.get(0).mStartPosition;
    mEndPosition = parts.get(parts.size() - 1).mEndPosition;

    final ArrayList<Long> longs = new ArrayList<>();
    for (AlignmentSection section : parts) {
      longs.add(section.mContig);
    }
    mContigs = longs;
  }
  /**
   * Represents the path through a graph that an alignment travels.
   * @param forward the <code>AlignmentChain</code> of contigs going forward
   * @param reverse the <code>AlignmentChain</code> of contigs going backward
   * @param graph the graph this alignment is against
   */
  public GraphAlignment(AlignmentChain forward, AlignmentChain reverse, Graph graph) {
    mScore = reverse.mScore + forward.mScore;
    mStartPosition = reverse.mSection.mStartPosition;
    mEndPosition = forward.mSection.mEndPosition;
    final ArrayList<Long> longs = new ArrayList<>();
    AlignmentChain current = reverse;
    while (current != null) {
      longs.add(current.mSection.mContig);
      current = current.mPrevious;
    }
    final int addPos = longs.size();
    current = forward;
    while (current.mPrevious != null) {
      longs.add(addPos, current.mSection.mContig);
      current = current.mPrevious;
    }
    mContigs = longs;
    mGraph = graph;

  }

  /**
   * Construct an alignment from through the provided contigs
   * @param startPosition start alignment position within the first contig
   * @param endPosition end alignment position within the last contig (inclusive)
   * @param contigs the contigs the alignment passes through
   * @param score the score of the alignment
   * @param graph the graph the alignment runs through
   */
  public GraphAlignment(int startPosition, int endPosition, List<Long> contigs, int score, Graph graph) {
    mStartPosition = startPosition;
    mEndPosition = endPosition;
    mContigs = contigs;
    mScore = score;
    mGraph = graph;
  }

  List<Long> contigs() {
    return mContigs;
  }
  long startContig() {
    return mContigs.get(0);
  }
  long endContig() {
    return mContigs.get(mContigs.size() - 1);
  }
  int startPosition() {
    return mStartPosition;
  }
  int endPosition() {
    return mEndPosition;
  }

  static boolean isPalindrome(long contigId, Graph graph) {
    final Contig contig = graph.contig(contigId);
    final Contig reverse = graph.contig(-contigId);
    for (int i = 0; i * 2 < contig.length(); ++i) {
      if (contig.nt(i) != reverse.nt(i)) {
        return false;
      }
    }
    return true;
  }

  @Override
  public boolean equals(Object o) {
    if (o == null || getClass() != o.getClass()) {
      return false;
    }

    final GraphAlignment that = (GraphAlignment) o;

    if (mScore != that.mScore) {
      return false;
    }
    if (startPosition() != that.startPosition()) {
      return false;
    }
    if (endPosition() != that.endPosition()) {
      return false;
    }

    if (contigs().size() != that.contigs().size()) {
      return false;
    }
    for (int i = 0; i < contigs().size(); ++i) {
      final long thisContig = contigs().get(i);
      final long thatContig = that.contigs().get(i);
      if (thisContig == thatContig) {
        continue;
      }
      if (thisContig == -thatContig && isPalindrome(thisContig, mGraph)) {
        continue;
      }
      return false;
    }
    return true;
  }

  @Override
  public int hashCode() {
    final int result =  Utils.pairHashContinuous(mScore, startPosition(), endPosition());
    long  contigSum = 0;
    for (long contig : contigs()) {
      contigSum += Math.abs(contig);
    }
    return Utils.pairHash(result, LongUtils.hashCode(contigSum));
  }
}
