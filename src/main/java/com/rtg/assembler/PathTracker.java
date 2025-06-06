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
    mPathCounts.merge(normalized, 1, (a, b) -> a + b);

  }

  List<Long> normalize(List<Long> path) {
    boolean forward = true;
    final int size = path.size();
    for (int i = 0; i < size; ++i) {
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
      for (int i = size - 1; i >= 0; --i) {
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
      for (int i = 0; i < o1.size() && i < o2.size(); ++i) {
        final int compare =  o1.get(i).compareTo(o2.get(i));
        if (compare != 0) {
          return compare;
        }
      }
      return Integer.compare(o1.size(), o2.size());

    }
  }

  static SortedMap<List<Long>, Integer> merge(List<PathTracker> trackers) {
    final SortedMap<List<Long>, Integer> sorted = new TreeMap<>(new PathComparator());
    for (PathTracker tracker : trackers) {
      for (Map.Entry<List<Long>, Integer> entry : tracker.mPathCounts.entrySet()) {
        final List<Long> path = entry.getKey();
        final Integer count = entry.getValue();
        sorted.merge(path, count, (a, b) -> b + a);
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
      graph.setPathAttribute(absId, GraphKmerAttribute.READ_COUNT, String.valueOf(count));
      totalReadsAdded += trackedReads;
    }
    return totalReadsAdded;
  }

}
