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
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.PathsIterator;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.implementation.PathArray;

/**
 */
@TestClass({"com.rtg.assembler.ConsensusTest", "com.rtg.assembler.MergeNodesTest"})
public class MergeNodes {
  final MutableGraph mGraph;
  final int mThreshold;
  final int mOverlap;

  /**
   * Construct a merger for the specified graph
   * @param graph the graph to merge in
   * @param threshold the number of reads required before a path is believed
   * @param overlap the number of bases adjacent contigs overlap by
   */
  public MergeNodes(MutableGraph graph, int threshold, int overlap) {
    mGraph = graph;
    mThreshold = threshold;
    mOverlap = overlap;
  }
  static final Set<String> INT_ATTRIBUTES = new LinkedHashSet<>();
  static {
    INT_ATTRIBUTES.add(GraphKmerAttribute.READ_COUNT);
    INT_ATTRIBUTES.add(GraphKmerAttribute.K_MER_FREQ);
  }

  static long mergeNodes(MutableGraph graph, int overlap, List<Long> mergedIds) {
    final long newId = graph.addContig(new ContigConcatenated(mergedIds, graph, overlap + 1));
    updatePaths(graph, newId, mergedIds);
    if (graph.contigAttributes().containsKey(Consensus.COMBINED)) {
      graph.setContigAttribute(newId, Consensus.COMBINED, ContigCollector.buildCombinedString(mergedIds, graph));
    }
    for (String attr : graph.contigAttributes().keySet()) {
      if (INT_ATTRIBUTES.contains(attr)) {
        int total = 0;
        for (long existingId : mergedIds) {
          final String stringValue = graph.contigAttribute(existingId, attr);
          if (stringValue != null) {
            total += Integer.parseInt(stringValue);
          }
        }
        if (total > 0) {
          graph.setContigAttribute(newId, attr, "" + total);
        }
      } else if (graph.contigAttribute(newId, attr) == null) {
        for (long existingId : mergedIds) {
          final String value = graph.contigAttribute(existingId, attr);
          if (value != null) {
            graph.setContigAttribute(newId, attr, value);
          }
        }
      }

    }

    return newId;
  }

  static void updatePaths(MutableGraph graph, long newId, List<Long> mergedIds) {
    final Set<Long> pathIds = new LinkedHashSet<>();
    for (long contig : mergedIds) {
      final PathsIterator it = graph.paths(contig);
      long pathId;
      while ((pathId = it.nextPathId()) != 0) {
        pathIds.add(pathId);
      }
    }
    for (long pathId : pathIds) {
      if (!graph.pathDeleted(pathId)) {
        final List<Long> existingPath = getPath(graph, pathId);
        final List<Long> newPath = newPath(newId, mergedIds, existingPath, graph);
        if (newPath.contains(newId)) {
          if (newPath.size() > 1) {
            final long matchingPath = GraphMap.findPath(newPath, graph, false);
            if (matchingPath == 0) {
              final long newPathId = graph.addPath(new PathArray(newPath));
              for (String key : graph.pathAttributes().keySet()) {
                final String value = graph.pathAttribute(pathId, key);
                if (value != null) {
                  graph.setPathAttribute(newPathId, key, value);
                }
              }
            } else if (graph.pathAttributes().containsKey(GraphKmerAttribute.READ_COUNT)) {
              final int newCount = GraphMap.existingReadCount(graph, matchingPath) + GraphMap.existingReadCount(graph, pathId);
              if (newCount > 0) {
                graph.setPathAttribute(matchingPath, GraphKmerAttribute.READ_COUNT, "" + newCount);

              }
            }
          }
          graph.deletePath(pathId);
        }
      }
    }
  }

  static List<Long> getPath(Graph graph, long pathId) {
    final ArrayList<Long> path = new ArrayList<>();
    for (int i = 0; i < graph.pathLength(pathId); ++i) {
      path.add(graph.pathContig(pathId, i));
    }
    return path;
  }

  static List<Long> newPath(long newNode, List<Long> mergedNodes, List<Long> existingPath, Graph graph) {
    final ArrayList<Long> newPath = new ArrayList<>();
    //Handle matches at start
    //Find longest overlap of mergedNodes with existingPath at the start of existingPath
    final int overlapSize = longestStartOverlap(mergedNodes, existingPath, graph);
    if (overlapSize > 0) {
      newPath.add(newNode);
    }
    for (int existingPos = overlapSize; existingPos < existingPath.size(); ++existingPos) {
      final boolean match = internalMatch(existingPos, mergedNodes, existingPath, graph);
      if (match) {
        newPath.add(newNode);
        existingPos += mergedNodes.size() - 1;
      } else {
        newPath.add(existingPath.get(existingPos));
      }
    }
    return newPath;
  }

  static boolean internalMatch(int start, List<Long> mergedNodes, List<Long> existingPath, Graph graph) {
    final int end = Math.min(mergedNodes.size(), existingPath.size() - start);
    final List<Long> mergedSub = mergedNodes.subList(0, end);
    final List<Long> existingSub = existingPath.subList(start, start + end);
    return palindromeSubMatch(graph, mergedSub, existingSub);
  }

  private static boolean palindromeSubMatch(Graph graph, List<Long> mergedSub, List<Long> existingSub) {
    for (int i = 0; i < mergedSub.size(); ++i) {
      final long merged = mergedSub.get(i);
      final long existing = existingSub.get(i);
      if (merged == existing) {
        continue;
      }
      if (merged == -existing && GraphAlignment.isPalindrome(merged, graph)) {
        continue;
      }
      return false;
    }
    return true;
  }


  static int longestStartOverlap(List<Long> mergedNodes, List<Long> existingPath, Graph graph) {
    for (int i = 0; i < mergedNodes.size(); ++i) {
      final int length = Math.min(mergedNodes.size() - i, existingPath.size());
      if (palindromeSubMatch(graph, mergedNodes.subList(i, i + length), existingPath.subList(0, length))) {
        return length;
      }
    }
    return 0;
  }

  void simplifyGraph() {
    for (long i = 1; i <= mGraph.numberContigs(); ++i) {
      if (!mGraph.contigDeleted(i)) {
        simplifyContig(i);
      }
      if (!mGraph.contigDeleted(-i)) {
        simplifyContig(-i);
      }
    }
  }

  private void simplifyContig(long current) {
    // We can only merge nodes when the first has a single successor.
    final Set<Long> successors = predecessors(mGraph, -current);
    if (isUniquePredecessor(successors) && !GraphAlignment.isPalindrome(current, mGraph)) {
      // use loop to get the single element due to absence of get method on sets
      // May also handle both orientations of a palindrome
      for (long successor : successors) {
        final List<Long> unambiguous = FilterPaths.unambiguousPath(current, -successor, mGraph, mThreshold);
        int i = unambiguous.size() - 1;
        while (i > 0 && !(hasUniquePredecessor(unambiguous.get(i)) && isUnambiguous(unambiguous, i))) {
          --i;
        }
        if (i > 0 && !GraphAlignment.isPalindrome(current, mGraph)) {
          final List<Long> mergeable = unambiguous.subList(0, i + 1);
//          System.err.println("Merging: " + mergeable);
          mergePath(mergeable);
        }
      }
    }
  }

  private boolean hasLinks(long current) {
    return predecessors(mGraph, current).size() != 0 || !predecessors(mGraph, -current).isEmpty();
  }
  boolean hasUniquePredecessor(long contig) {
    final Set<Long> predecessors = predecessors(mGraph, contig);
    return isUniquePredecessor(predecessors);
  }
  private boolean isUniquePredecessor(Set<Long> predecessors) {
    if (predecessors.size() <= 1) {
      return true;
    } else if (predecessors.size() == 2) {
      final Set<Long> forward = new LinkedHashSet<>(predecessors.size());
      for (long l : predecessors) {
        forward.add(Math.abs(l));
      }
      if (forward.size() != 1) {
        return false;
      }
      for (long l : forward) {
        if (!GraphAlignment.isPalindrome(l, mGraph)) {
          return false;

        }
      }
      return true;
    } else {
      return false;
    }
  }

  private void mergePath(List<Long> contigs) {
    final ArrayList<Boolean> hadLinks = new ArrayList<>();
    for (long contig : contigs) {
      hadLinks.add(hasLinks(contig));
    }
    mergeNodes(mGraph, mOverlap, contigs);
    for (int i = 0; i < hadLinks.size(); ++i) {
      final long current = contigs.get(i);
      if (hadLinks.get(i) && !hasLinks(current)) {
        mGraph.deleteContig(current);
      }
    }
  }

  private boolean isUnambiguous(List<Long> contigs, int end) {
    final List<Long> reverse = FilterPaths.unambiguousPath(-contigs.get(end), -contigs.get(end - 1), mGraph, mThreshold);
    if (reverse.size() <= end) {
      return false;
    }
    if (reverse.get(end) != -contigs.get(0)) {
      return false;
    }
    return true;

  }

  static Set<Long> predecessors(MutableGraph graph, long i) {
    final Set<Long> predecessors = new LinkedHashSet<>();
    final PathsIterator iterator = graph.paths(i);
    long pathId;
    while ((pathId = iterator.nextPathId()) != 0) {
      if (iterator.contigIndex() > 0) {
        predecessors.add(graph.pathContig(pathId, iterator.contigIndex() - 1));
      }
    }
    return predecessors;
  }
}
