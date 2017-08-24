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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.rtg.assembler.graph.Graph;
import com.rtg.util.LongUtils;
import com.rtg.util.Utils;

/**
 */
public class PairJoiner {

  private static final int MAX_PATHS_PER_START_POSITION = 1000; //Integer.parseInt(System.getProperty("rtg.assembler.maxpaths", "1000"));
  static final Set<GraphAlignment> TOO_MANY_PATHS = Collections.emptySet();

  final Graph mGraph;
  final GraphTraversions mTraversions;
  final int mOverlap;

  /**
   * Finds paths between potential alignment positions in the graph
   * @param graph the graph to explore
   * @param overlap number of bases which overlap in adjacent contigs
   */
  public PairJoiner(Graph graph, int overlap) {
    this(graph, overlap, new GraphTraversions(graph));
  }
  /**
   * Finds paths between potential alignment positions in the graph
   * @param graph the graph to explore
   * @param overlap number of bases which overlap in adjacent contigs
   * @param traverse the links between contigs in the graph that you will consider
   */
  public PairJoiner(Graph graph, int overlap, GraphTraversions traverse) {
    mGraph = graph;
    mTraversions = traverse;
    mOverlap = overlap;
  }

  Set<GraphAlignment> paired(List<Set<GraphAlignment>> fragmentAlignments, int minInsert, int maxInsert) {
    if (fragmentAlignments.size() < 1) {
      return Collections.emptySet();
    }
    Set<GraphAlignment> fromSet = fragmentAlignments.get(0);
    for (int i = 1; i < fragmentAlignments.size(); ++i) {
      final Set<GraphAlignment> toSet = fragmentAlignments.get(i);
      final Map<Long, Set<GraphAlignment>> matchMap = destinationMap(toSet, mGraph);
      final Set<GraphAlignment> joined = new HashSet<>();
      for (GraphAlignment start : fromSet) {
        final Set<GraphAlignment> alignments = joinAlignments(start, matchMap, minInsert, maxInsert);
        if (alignments == TOO_MANY_PATHS) {
          return alignments; // We are bailing out, so stop searching
        }
        joined.addAll(alignments);
      }
      fromSet = joined;
    }
    return fromSet;
  }

  static Map<Long, Set<GraphAlignment>> destinationMap(Set<GraphAlignment> toSet, Graph graph) {
    final Map<Long, Set<GraphAlignment>> matchMap = new HashMap<>();
    for (GraphAlignment alignment : toSet) {
      final long key = -alignment.endContig();
      if (matchMap.containsKey(key)) {
        matchMap.get(key).add(alignment);
      } else {
        final Set<GraphAlignment> list = new HashSet<>();
        list.add(alignment);
        matchMap.put(key, list);
      }
    }
    final List<Long> palindromes = new ArrayList<>();
    // add any alignments ending in a palindrome to the - strand
    for (long contigId : matchMap.keySet()) {
      if (GraphAlignment.isPalindrome(contigId, graph)) {
        palindromes.add(contigId);
      }
    }
    for (long contigId : palindromes) {
      final Set<GraphAlignment> alignments = matchMap.get(contigId);
      if (matchMap.containsKey(-contigId)) {
        matchMap.get(-contigId).addAll(alignments);
      } else {
        matchMap.put(-contigId, alignments);
      }
    }
    return matchMap;
  }

  static class PairChain {
    PairChain mPrevious;
    int mInsertLength;
    long mContigId;

    PairChain(PairChain previous, int length, long contigId) {
      mPrevious = previous;
      mInsertLength = length;
      mContigId = contigId;
    }

    @Override
    public String toString() {
      return "PairChain{" + "mPrevious=" + mPrevious + ", mInsertLength=" + mInsertLength + ", mContigId=" + mContigId + '}';
    }

    @Override
    public boolean equals(Object o) {
      if (o == null || getClass() != o.getClass()) {
        return false;
      }
      final PairChain that = (PairChain) o;
      if (mInsertLength != that.mInsertLength) {
        return false;
      }
      if (mContigId != that.mContigId) {
        return false;
      }
      return mPrevious == that.mPrevious;
    }


    @Override
    public int hashCode() {
      // Hashcode independent of previous to avoid long recursive calculation
      // Need this hashcode because we put these in a set
      return Utils.pairHash(mInsertLength, LongUtils.hashCode(mContigId));
    }
  }

  private boolean checkPathAgrees(final GraphAlignment start, final GraphAlignment destination, final int startPos) {
    for (int j = startPos, k = destination.contigs().size() - 1; k >= 0 && j < start.contigs().size(); ++j, --k) {
      final long forwardContig = start.contigs().get(j);
      final long reverseContig = destination.contigs().get(k);
      // Check match allowing for differing direction and presence of palindromes
      if (forwardContig != -reverseContig
            && (forwardContig != reverseContig || !GraphAlignment.isPalindrome(reverseContig, destination.mGraph))) {
        return false;
      }
    }
    return true;
  }

  Set<GraphAlignment> joinAlignments(GraphAlignment start, Map<Long, Set<GraphAlignment>> toSet, int minInsert, int maxInsert) {
    final Set<GraphAlignment> results = new HashSet<>();
    final PairChain initial = constructInitial(start);
    Set<PairChain> paths = new HashSet<>();
    paths.add(initial);
    // Resolve cross contig sloptigs
    for (int i = 0; i < start.contigs().size() - 1; ++i) {
      final long contig = start.contigs().get(i);
      final Set<GraphAlignment> destinations = toSet.get(contig);
      if (destinations != null) {
        for (GraphAlignment destination : destinations) {
          if (checkPathAgrees(start, destination, i)) {
            final int insertSize = slopSize(start, i, destination);
            if (insertSize < maxInsert && insertSize >= minInsert) {
              final GraphAlignment alignment = buildSloptigAlignment(start, i, destination);
              results.add(alignment);
            }
          }
        }
      }
    }
    // Non sloptigs
    while (paths.size() > 0) {
      if (paths.size() > MAX_PATHS_PER_START_POSITION) {
        return TOO_MANY_PATHS;
      }
      final Set<PairChain> next = new HashSet<>();
      for (PairChain chain : paths) {
        final Set<GraphAlignment> destinations = toSet.get(chain.mContigId);
        if (destinations != null) {
          for (GraphAlignment alignment : destinations) {
            final int insertSize = insert(chain, alignment);
            if (insertSize < maxInsert && insertSize >= minInsert) {
              results.add(buildPairAlignment(start.startPosition(), start.mScore, chain, alignment));
            }
          }
        }
        final int newInserLength = chain.mInsertLength + mGraph.contigLength(chain.mContigId);
        if (newInserLength - mOverlap >= maxInsert) {
          continue;
        }
        final Traversion traversion = mTraversions.get(chain.mContigId);
        for (long contigId : traversion.next()) {
          next.add(new PairChain(chain, newInserLength - mOverlap, contigId));
        }
      }
      paths = next;
    }

    return results;
  }

  private int slopSize(GraphAlignment start, int overlapContig, GraphAlignment destination) {
    int insert = mGraph.contigLength(destination.endContig()) - destination.endPosition() - 1;
    for (int i = overlapContig; i < start.contigs().size() - 1; ++i) {
      insert -= mGraph.contigLength(start.contigs().get(i)) - mOverlap;
    }
    insert -= start.endPosition();
    return insert;
  }

  private GraphAlignment buildSloptigAlignment(GraphAlignment startAlignment , int overlapContig, GraphAlignment endAlignment) {
    final int endPosition = mGraph.contigLength(endAlignment.startContig()) - endAlignment.startPosition() - 1;
    final List<Long> contigs = new ArrayList<>();
    contigs.addAll(startAlignment.contigs().subList(0, overlapContig));
    for (long contig : endAlignment.contigs()) {
      contigs.add(overlapContig, -contig);
    }
    return new GraphAlignment(startAlignment.startPosition(), endPosition, contigs, startAlignment.mScore + endAlignment.mScore, mGraph);
  }
  private GraphAlignment buildPairAlignment(int startPosition, int startScore, PairChain chain, GraphAlignment alignment) {
    assert chain != null;
    final List<Long> contigs = new ArrayList<>();
    PairChain current = chain;
    while (current != null) {
      contigs.add(0, current.mContigId);
      current = current.mPrevious;
    }
    final int insertPos = contigs.size();
    for (long contigId : alignment.contigs()) {
      contigs.add(insertPos, -contigId);
    }
    contigs.remove(insertPos);
    final int endPosition = mGraph.contigLength(alignment.startContig()) - alignment.startPosition() - 1;
    return new GraphAlignment(startPosition, endPosition, contigs, startScore + alignment.mScore, mGraph);
  }

  private int insert(PairChain chain, GraphAlignment alignment) {
    return chain.mInsertLength + mGraph.contigLength(chain.mContigId) - alignment.endPosition() - 1;
  }

  private PairChain constructInitial(GraphAlignment start) {
    PairChain previous = null;
    final int score = -start.endPosition() - 1;
    for (long contigId : start.contigs()) {
      previous = new PairChain(previous, score, contigId);
    }
    if (previous == null) {
      throw new IllegalArgumentException();
    }

    return previous;
  }
}
