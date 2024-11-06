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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.rtg.util.Pair;

/**
 */
public class ConstraintCache {
  Map<Long, List<ConstraintCollector>> mIndex = new HashMap<>();
  Map<Pair<Long, Long>, ConstraintCollector> mCache = new HashMap<>();

  List<ConstraintCollector> find(long contigId) {
    return mIndex.get(contigId);
  }
  ConstraintCollector find(long contigId, long mate) {
    final ConstraintCollector constraintCollector = mCache.get(normalizedPair(contigId, mate));
    if (constraintCollector == null) {
      return constraintCollector;
    }
    if (contigId < mate) {
      return constraintCollector;
    } else {
      return reverse(constraintCollector);
    }
  }

  static ConstraintCollector reverse(ConstraintCollector forward) {
    final ConstraintCollector reverse = new ConstraintCollector(forward.mContigB, forward.mContigA);
    reverse.mConstraint.addAll(forward.mConstraint);
    return reverse;
  }

  static ConstraintCache combineCaches(List<ConstraintCache> caches) {
    final ConstraintCache combined = new ConstraintCache();
    for (ConstraintCache cache : caches) {
      for (Map.Entry<Pair<Long, Long>, ConstraintCollector> entry : cache.mCache.entrySet()) {
        final Pair<Long, Long> key = entry.getKey();
        final ConstraintCollector value = entry.getValue();
        for (int i = 0; i < value.mConstraint.size(); ++i) {
          combined.addConstraint(key.getA(), key.getB(), 0, 0, value.mConstraint.get(i));
        }
      }
    }

    return combined;
  }

  Pair<Long, Long> normalizedPair(long a, long b) {
    if (a < b) {
      return Pair.create(a, b);
    } else {
      return Pair.create(b, a);
    }
  }
  void addConstraint(long contig1, long contig2, int distanceFromEnd1, int distanceFromEnd2, int expected) {
    final Pair<Long, Long> normal = normalizedPair(contig1, contig2);
    ConstraintCollector collector = mCache.get(normal);
    if (collector == null) {
      collector = new ConstraintCollector(normal.getA(), normal.getB());
      mCache.put(normal, collector);
      addToIndex(collector);
    }
    if (contig1 < contig2) {
      collector.increment(distanceFromEnd1, distanceFromEnd2, expected);
    } else {
      collector.increment(distanceFromEnd2, distanceFromEnd1, expected);
    }
  }
  void addToIndex(ConstraintCollector constraint) {
    addSingle(constraint.mContigA, constraint);
    addSingle(constraint.mContigB, constraint);
  }
  void addSingle(long id, ConstraintCollector constraint) {
    mIndex.computeIfAbsent(id, k -> new ArrayList<>()).add(constraint);
  }
}
