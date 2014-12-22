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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.implementation.PathArray;

/**
 */
public class PathTracker {
  final PalindromeTracker mPalindromes;
  final Map<List<Long>, Integer> mPathCounts = new HashMap<>();

  /**
   * Construct a path tracker taking the supplied palindromes into account
   * @param palindromes list of palindromes in the graph
   */
  public PathTracker(PalindromeTracker palindromes) {
    mPalindromes = palindromes;
  }

  static long findOrAddPath(List<Long> contigs, MutableGraph graph) {
    final long existing = GraphMap.findPath(contigs, graph, true);
    final long result;
    if (existing == 0) {
      result = graph.addPath(new PathArray(contigs));
    } else {
      result = existing;
    }
    return result;
  }

  void increment(List<Long> path) {
    final List<Long> normalized = normalize(path);
    final Integer existing = mPathCounts.get(normalized);
    if (existing == null) {
      mPathCounts.put(normalized, 1);
    } else {
      mPathCounts.put(normalized, existing + 1);
    }

  }

  List<Long> normalize(List<Long> path) {
    boolean forward = true;
    final int size = path.size();
    for (int i = 0; i < size; i++) {
      long f = path.get(i);
      f = f < 0 && mPalindromes.isPalindrome(f) ? Math.abs(f) : f;
      long r = -path.get(size - i - 1);
      r = r < 0 && mPalindromes.isPalindrome(r) ? Math.abs(r) : r;
      if (f < r) {
        forward = true;
        break;
      } else if (f > r) {
        forward = false;
        break;
      }
    }
    final List<Long> result = new ArrayList<>();
    if (forward) {
      for (Long aPath : path) {
        long f = aPath;
        f = f < 0 && mPalindromes.isPalindrome(f) ? Math.abs(f) : f;
        result.add(f);
      }
    } else {
      for (int i = size - 1; i >= 0; i--) {
        long f = -path.get(i);
        f = f < 0 && mPalindromes.isPalindrome(f) ? Math.abs(f) : f;
        result.add(f);
      }
    }
    return result;
  }
  static class PathComparator implements Comparator<List<Long>>, Serializable {

    @Override
    public int compare(List<Long> o1, List<Long> o2) {
      for (int i = 0; i < o1.size() && i < o2.size(); i++) {
        final int compare =  o1.get(i).compareTo(o2.get(i));
        if (compare != 0) {
          return compare;
        }
      }
      return Integer.valueOf(o1.size()).compareTo(o2.size());

    }
  }

  static SortedMap<List<Long>, Integer> merge(List<PathTracker> trackers) {
    final SortedMap<List<Long>, Integer> sorted = new TreeMap<>(new PathComparator());
    for (PathTracker tracker : trackers) {
      for (Map.Entry<List<Long>, Integer> entry : tracker.mPathCounts.entrySet()) {
        final List<Long> path = entry.getKey();
        final Integer count = entry.getValue();
        final Integer existing = sorted.get(path);
        if (existing == null) {
          sorted.put(path, count);
        } else {
          sorted.put(path, count + existing);
        }
      }
    }
    return sorted;
  }
  static int apply(SortedMap<List<Long>, Integer> merged, MutableGraph graph) {
    int totalReadsAdded = 0;
    for (Map.Entry<List<Long>, Integer> entry: merged.entrySet()) {
      final long pathId = findOrAddPath(entry.getKey(), graph);
      final long absId = Math.abs(pathId);
      final String readCount = graph.pathAttribute(absId, GraphKmerAttribute.READ_COUNT);
      int count;
      if (readCount == null) {
        count = 0;
      } else {
        count = Integer.parseInt(readCount);
      }
      final int trackedReads = entry.getValue();
      count += trackedReads;
      graph.setPathAttribute(absId, GraphKmerAttribute.READ_COUNT, "" + count);
      totalReadsAdded += trackedReads;
    }
    return totalReadsAdded;
  }

}
