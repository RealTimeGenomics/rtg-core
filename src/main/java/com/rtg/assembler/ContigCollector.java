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

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.implementation.IteratorLocal;
import com.rtg.util.array.intindex.IntChunks;

/**
 */
public class ContigCollector {
  final long mTipThreshold;
  final IntChunks mStartTips;
  final IntChunks mEndTips;
  final GraphKmerAttribute mGraph;
  final IteratorLocal mPathIterator;
  private final int mKmerSize;

  /**
   * Class that collapses un-ranching strings of contigs into a single contig. It ignores branches to tip when there is
   * a non-tip alternative
   * @param tipThreshold maximum tip value below which branches will be ignored
   * @param kmerSize length of a kmer. Adjacent contigs overlap by 1 less than this value.
   * @param startTips start tip values for the contigs in <code>graph</code>
   * @param endTips end tip values for the contigs in <code>graph</code>
   * @param graph graph implementation containing the contigs
   */
  public ContigCollector(long tipThreshold, int kmerSize, IntChunks startTips, IntChunks endTips, GraphKmerAttribute graph) {
    mTipThreshold = tipThreshold;
    mKmerSize = kmerSize;
    mStartTips = startTips;
    mEndTips = endTips;
    mGraph = graph;
    mPathIterator = new IteratorLocal(graph);
    if (mStartTips == null ^ mEndTips == null) {
      throw new IllegalArgumentException();
    }
    assert mStartTips == null || mStartTips.length() == mEndTips.length();
  }

  boolean tipValueOk(long id) {
    final long i = Math.abs(id);
    if (mStartTips == null || i >= mStartTips.length()) {
      return true;
    }
    final int startTip = (int) mStartTips.get(i);
    final int endTip = (int) mEndTips.get(i);
    return !(startTip != 0 && startTip < mTipThreshold) && !(endTip != 0 && endTip < mTipThreshold);
  }

  void collapse() {
    for (long i = 1; i <= mGraph.numberContigs(); ++i) {
      if (!tipValueOk(i) || mGraph.contigDeleted(i)) {
        continue;
      }
      final List<Long> walkedIds = new ArrayList<>();
      walkedIds.add(i);
      walk(walkedIds, i, true);
      walk(walkedIds, i, false);
//      System.err.println(walkedIds);
      if (walkedIds.size() > 1) {
        final long newId = MergeNodes.mergeNodes(mGraph, mKmerSize - 1, walkedIds);
        if (mGraph.contigAttributes().containsKey(Consensus.COMBINED)) {
          final String combined = buildCombinedString(walkedIds, mGraph);
          mGraph.setContigAttribute(newId, Consensus.COMBINED, combined);
        }
        if (mStartTips != null) {
          newTipValues(walkedIds, newId);
        }
        int kmerCount = 0;
        for (long id : walkedIds) {
          kmerCount += mGraph.kmerFreq(id);
          mGraph.deleteContig(id);
        }
        mGraph.setKmerFreq(newId, kmerCount);
      }
    }
  }

  static String buildCombinedString(List<Long> walkedIds, Graph graph) {
    final StringBuilder sb = new StringBuilder();
    String join = "";
    for (long id : walkedIds) {
      sb.append(join);
      final String combined = graph.contigAttribute(id, Consensus.COMBINED);
      sb.append(id);
      if (combined != null) {
        sb.append(":");
        sb.append("(").append(combined).append(")");
      }
      join = "/";
    }
    return sb.toString();
  }

  private void newTipValues(List<Long> walkedIds, long newId) {
    // New tip value should be the largest tip value of the combined set or 0 if any of the contigs have that.
    // We've got to account for RC links
    int startTip = 0;
    int endTip = 0;
    // track if something is not a tip, if these end up true then the whole tip value will be 0
    boolean zeroStart = false;
    boolean zeroEnd = false;
    for (long id : walkedIds) {
      final int endTipVal = Math.abs(id) >= mEndTips.length() ? 0 : mEndTips.getInt(Math.abs(id));
      final int startTipVal = Math.abs(id) >= mStartTips.length() ? 0 : mStartTips.getInt(Math.abs(id));
      final int start = id < 0 ? endTipVal : startTipVal;
      final int end = id < 0 ? startTipVal : endTipVal;
      // included forward so swap
      zeroStart |= start == 0;
      zeroEnd |= end == 0;
      startTip = Math.max(startTip, start);
      endTip = Math.max(endTip, end);
    }
    mStartTips.extendTo(newId + 1);
    mStartTips.set(newId, zeroStart ? 0 : startTip);
    mEndTips.extendTo(newId + 1);
    mEndTips.set(newId, zeroEnd ? 0 : endTip);
  }

  void walk(List<Long> walked, long start, boolean direction) {
    long current = start;
    while ((current = monogamousLink(current, direction)) != 0) {
      if (walked.contains(current)) {
        break;
      }
      if (direction) {
        walked.add(current);
      } else {
        walked.add(0, current);
      }
    }

  }

  /**
   * @param id current <code>contigId</code>
   * @param direction true if forward, false if reverse
   * @return the contig id of the single next contig in <code>direction</code> provided it only has <code>id</code> as
   *         a previous contig
   */
  long monogamousLink(long id, boolean direction) {
    final long next = uniqueNext(id, direction);
    if (next == 0) {
      return 0;
    }
    final long reverse = uniqueNext(next, !direction);
    return reverse == id ? next : 0;
  }

  /**
   * @param id current contig id
   * @param forward true if forward, false if reverse
   * @return the contig id of the single next contig in <code>forward</code>
   */
  long uniqueNext(long id, boolean forward) {
    mPathIterator.set(id);
    long next = 0;
    long pathId;
    boolean seenTips = false;
    while ((pathId = mPathIterator.nextPathId()) != 0) {
      final int index = mPathIterator.contigIndex();
      final int nextIndex = index + (forward ? 1 : -1);
      if (nextIndex > mGraph.pathLength(pathId) - 1 || nextIndex < 0) {
        continue;
      }
      final long linkedContig = mGraph.pathContig(pathId, nextIndex);
      if (next != 0 && next != linkedContig) {
        if (tipValueOk(next) && tipValueOk(linkedContig)) {
          // multiple non-tip paths
          return 0;
        } else if (!tipValueOk(next) && tipValueOk(linkedContig)) {
          // first non-tip contig
          next = linkedContig;
        } else if (!tipValueOk(next) && !tipValueOk(linkedContig)) {
          // Encountered multiple tip links, need to wait for a non-tip
          seenTips = true;
          next = 0;
        }
      } else {
        if (tipValueOk(linkedContig) || !seenTips) {
          next = linkedContig;
        }
      }
    }
    return next;
  }

}
