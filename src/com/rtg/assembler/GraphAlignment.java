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
