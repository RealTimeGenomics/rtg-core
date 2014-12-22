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
    for (int i = 0; i < forward.mConstraint.size(); i++) {
      reverse.mConstraint.add(forward.mConstraint.get(i));
    }
    return reverse;
  }

  static ConstraintCache combineCaches(List<ConstraintCache> caches) {
    final ConstraintCache combined = new ConstraintCache();
    for (ConstraintCache cache : caches) {
      for (Map.Entry<Pair<Long, Long>, ConstraintCollector> entry : cache.mCache.entrySet()) {
        final Pair<Long, Long> key = entry.getKey();
        final ConstraintCollector value = entry.getValue();
        for (int i = 0; i < value.mConstraint.size(); i++) {
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
    List<ConstraintCollector> list =  mIndex.get(id);
    if (list == null) {
      list = new ArrayList<>();
      mIndex.put(id, list);
    }
    list.add(constraint);
  }
}
