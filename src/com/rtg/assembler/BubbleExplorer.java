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
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.rtg.assembler.graph.PathsIterator;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.implementation.IteratorLocal;
import com.rtg.assembler.graph.implementation.PathArray;
import com.rtg.util.Histogram;

/**
 *         Date: 8/05/12
 *         Time: 2:08 PM
 */
public class BubbleExplorer {
  private final GraphKmerAttribute mGraph;
  private static final int COMPLEXITY_THRESHOLD = 100;
  private final double mMergeThreshold;
  private final int mMaxBubbleSize;
  private final int mMaxBubbleDifference;
  private final int mKmerSize;
  private final IteratorLocal mPathsIterator;

  static class PendingNode {
    long mNodeId;
    long mPreviousNodeId;
    int mCumulativeCoverage;
    List<Long> mLinksIn;
    int mCountLinksInExplored;
    List<Long> mLinksOut;
    int mMinLength;
    int mMaxLength;
    int mLength;
    int mNodeCoverage;


    PendingNode(long nodeId, long previousNodeId, int cumulativeCoverage, GraphKmerAttribute graph, IteratorLocal pathsIterator) {
      mNodeId = nodeId;
      mPreviousNodeId = previousNodeId;
      mCumulativeCoverage = cumulativeCoverage;
      mLength = graph.contigLength(nodeId);
      mMinLength = mLength;
      mMaxLength = mLength;
      mCountLinksInExplored = 1;
      pathsIterator.set(nodeId);
      mLinksIn = links(false, graph, pathsIterator);
      pathsIterator.set(nodeId);
      mLinksOut = links(true, graph, pathsIterator);
      mNodeCoverage = graph.kmerFreq(nodeId);
    }


    boolean allLinksInFollowed() {
      return mLinksIn.size() == mCountLinksInExplored;
    }

    @Override
    public String toString() {
      return "PendingNode{"
          + "mNodeId=" + mNodeId
          + " mMinLength=" + mMinLength
          + " mMaxLength=" + mMaxLength
          + '}';
    }
  }

  static List<Long> links(boolean end, GraphKmerAttribute graph, PathsIterator pathsIterator) {
    final List<Long> links = new ArrayList<>();
    long pathId;
    while ((pathId = pathsIterator.nextPathId()) != 0) {
      final int index = pathsIterator.contigIndex();
      final long linkedContig = graph.pathContig(pathId, 1 - index);
      if (end && index == 0) {
        links.add(linkedContig);
      } else if (!end && index == 1) {
        links.add(linkedContig);
      }
    }
    return links;
  }

  BubbleExplorer(GraphKmerAttribute graph, int kmerSize, int maxBubbleSize, int maxBubbleDifference) {
    this(graph, kmerSize, maxBubbleSize, maxBubbleDifference, 0.0);
  }

  BubbleExplorer(GraphKmerAttribute graph, int kmerSize, int maxBubbleSize, int maxBubbleDifference, double mergeThreshold) {
    mGraph = graph;
    mMaxBubbleSize = maxBubbleSize;
    mMaxBubbleDifference = maxBubbleDifference;
    mKmerSize = kmerSize;
    mPathsIterator = new IteratorLocal(graph);
    mMergeThreshold = mergeThreshold;
  }
  Histogram popBubbles() {
    final Histogram h = new Histogram();
    final List<Long> mergeCandidates = new ArrayList<>();
    for (long i = 1; i <= mGraph.numberContigs(); ++i) {
      mPathsIterator.set(i);
      if (links(true, mGraph, mPathsIterator).size() > 1) {
        mergeCandidates.add(i);
      }
      mPathsIterator.set(i);
      if (links(false, mGraph, mPathsIterator).size() > 1) {
        mergeCandidates.add(-i);
      }
    }
    final Set<Long> removed = new HashSet<>();
    for (int i = 0; i < mergeCandidates.size(); ++i) {
      final long nodeId = mergeCandidates.get(i);
      if (removed.contains(nodeId)) {
        continue;
      }
      boolean success = true;
      while (success) {
        final Bubble b = explore(nodeId);
        if (b.mResult == BubbleResult.FOUND) {
          int kmerCount = 0;
          for (long mergedId : b.mMergedNodes) {
            kmerCount += mGraph.kmerFreq(mergedId);
          }
//          System.err.println(b.mKmerRatio + " " + mMergeThreshold);
          if (b.mKmerRatio >= mMergeThreshold) {
            h.increment(b.mMergedNodes.size());
            final ContigConcatenated contigConcatenated = b.mergedNode(mGraph, mKmerSize, mPathsIterator);
            final long newId = mGraph.addContig(contigConcatenated);
            for (long link : contigConcatenated.getStartLinks()) {
              mGraph.addPath(new PathArray(link, newId));
            }
            for (long link : contigConcatenated.getEndLinks()) {
              mGraph.addPath(new PathArray(newId, link));
            }
            for (long mergedId : b.mMergedNodes) {
              mGraph.deleteContig(mergedId);
              removed.add(mergedId);
              removed.add(-mergedId);
            }
            mGraph.setKmerFreq(newId, kmerCount);
            mergeCandidates.add(newId);
            mergeCandidates.add(-newId);
          } else {
            success = false;
          }
        } else {
          success = false;
        }
      }
    }
    return h;
  }

  Bubble explore(long node) {
    final PendingNode startNode = new PendingNode(node, 0, mGraph.kmerFreq(node), mGraph, mPathsIterator);


    final Map<Long, PendingNode> pendingNodes = new HashMap<>();
    final Map<Long, PendingNode> exploredNodes = new HashMap<>();
    exploredNodes.put(startNode.mNodeId, startNode);
    for (Long next : startNode.mLinksOut) {
      final PendingNode primeNode = new PendingNode(next, startNode.mNodeId, startNode.mCumulativeCoverage + mGraph.kmerFreq(next), mGraph, mPathsIterator);
      pendingNodes.put(next, primeNode);
    }
    // Loop until paths all converge
    bubbleLoop: while (pendingNodes.size() > 1) {
      if (pendingNodes.size() > COMPLEXITY_THRESHOLD) {
        break;
      }

      for (PendingNode pending : pendingNodes.values()) {
        // no nodes all linksIn explored and no linksOut
        if (pending.mLinksOut.isEmpty() && pending.allLinksInFollowed()) {
          return IMPOSSIBLE_BUBBLE;
        }
      }

      //get first expandable node
      PendingNode expandMe = null;
      for (PendingNode pending : pendingNodes.values()) {
        if (pending.allLinksInFollowed()) {
          expandMe = pending;
          break;
        }
      }
      if (expandMe == null) {
        return IMPOSSIBLE_BUBBLE;
      }
      // no nodes with maxLength > maxBubbleSize
      if (expandMe.mMaxLength > mMaxBubbleSize) {
        break;
      }
      for (Long next : expandMe.mLinksOut) {

        if (exploredNodes.containsKey(next) || exploredNodes.containsKey(-next)) {
          // Loop?
          return IMPOSSIBLE_BUBBLE;
        }
        if (pendingNodes.containsKey(next)) {
          final PendingNode existing = pendingNodes.get(next);
          final int nodeLength = existing.mLength - mKmerSize + 1;
          existing.mMaxLength = Math.max(existing.mMaxLength, expandMe.mMaxLength + nodeLength);
          existing.mMinLength = Math.min(existing.mMinLength, expandMe.mMinLength + nodeLength);
          if (existing.mCumulativeCoverage < expandMe.mCumulativeCoverage + existing.mNodeCoverage) {
            existing.mPreviousNodeId = expandMe.mNodeId;
            existing.mCumulativeCoverage = expandMe.mCumulativeCoverage + existing.mNodeCoverage;
          }
          existing.mCountLinksInExplored++;
          // no nodes with maxLength - minLength > maxBubbleDifference
          if (existing.mMaxLength - existing.mMinLength > mMaxBubbleDifference) {
            break bubbleLoop;
          }
        } else {
          final PendingNode created = new PendingNode(next, expandMe.mNodeId, expandMe.mCumulativeCoverage + mGraph.kmerFreq(next), mGraph, mPathsIterator);
          created.mMaxLength += expandMe.mMaxLength - mKmerSize + 1;
          created.mMinLength += expandMe.mMinLength - mKmerSize + 1;
          pendingNodes.put(next, created);
        }
      }
      pendingNodes.remove(expandMe.mNodeId);
      exploredNodes.put(expandMe.mNodeId, expandMe);
    }

    if (pendingNodes.size() == 1) {
      final PendingNode finalNode = pendingNodes.values().iterator().next();
      PendingNode current = finalNode;
      if (current.allLinksInFollowed()) {
        int bestKmerFreq = 0;
        final List<Long> bestPath = new ArrayList<>();
        bestPath.add(current.mNodeId);
        while (current.mPreviousNodeId != 0) {
          current = exploredNodes.get(current.mPreviousNodeId);
          bestKmerFreq += mGraph.kmerFreq(current.mNodeId);
          bestPath.add(0, current.mNodeId);
        }
        // excluding first and last kmer freqs because they're in all paths and will outweigh the branches
        bestKmerFreq -= mGraph.kmerFreq(current.mNodeId);
        final List<Long> mergedNodes = new ArrayList<>();
        mergedNodes.add(finalNode.mNodeId);
        mergedNodes.addAll(exploredNodes.keySet());
        int allKmerFreq = 0;
//        System.err.println(mergedNodes);
        for (long merged : mergedNodes) {
//          System.err.print(mGraph.kmerFreq(merged) + ", ");
          if (merged != finalNode.mNodeId && merged != current.mNodeId) {
            allKmerFreq += mGraph.kmerFreq(merged);
          }
        }
//        System.err.println();
//        System.err.println("bestKmerFreq:" + bestKmerFreq + " allKmerFreq:" + allKmerFreq);
        final double ratio = allKmerFreq == 0 ? 1.0 : (double) bestKmerFreq / (double) allKmerFreq;
        return new Bubble(bestPath, mergedNodes, BubbleResult.FOUND, ratio);
      }
    }
    // failed to find a path
    return NOT_FOUND_BUBBLE;
  }
  static final Bubble IMPOSSIBLE_BUBBLE = new Bubble(null, null, BubbleResult.IMPOSSIBLE, 0.0);
  static final Bubble NOT_FOUND_BUBBLE = new Bubble(null, null, BubbleResult.NOT_FOUND, 0.0);
  enum BubbleResult {
    FOUND
    , NOT_FOUND
    , IMPOSSIBLE
  }
  static class Bubble {
    List<Long> mBestPath;
    Collection<Long> mMergedNodes;
    final BubbleResult mResult;
    final double mKmerRatio;

    Bubble(List<Long> bestPath, Collection<Long> mergedNodes, BubbleResult result, double bestRatio) {
      mBestPath = bestPath;
      mMergedNodes = mergedNodes;
      mResult = result;
      mKmerRatio = bestRatio;
    }
    ContigConcatenated mergedNode(GraphKmerAttribute graph, int kmerSize, IteratorLocal pathsIterator) {
      final ContigConcatenated concat = new ContigConcatenated(mBestPath, graph, kmerSize);
      final long start = mBestPath.get(0);
      final long end = mBestPath.get(mBestPath.size() - 1);
      pathsIterator.set(end);
      concat.setEndLinks(links(true, graph, pathsIterator));
      pathsIterator.set(start);
      concat.setStartLinks(links(false, graph, pathsIterator));

      return concat;
    }

    @Override
    public String toString() {
      return "Bubble{"
          + "mBestPath=" + mBestPath
          + ", mMergedNodes=" + mMergedNodes
          + ", mResult=" + mResult
          + '}';
    }
  }

}
