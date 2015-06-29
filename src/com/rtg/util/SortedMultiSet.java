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

package com.rtg.util;

import java.util.Map.Entry;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * An ordered multi-set.
 * @param <E> type of elements in the multi-set.
 */
public class SortedMultiSet<E extends Comparable<E>> extends MultiSet<E> {

  /**
   * Constructor.
   */
  public SortedMultiSet() {
    super(new TreeMap<E, Counter>());
  }

  /**
   * @return all the entries in the multi-set in ascending order of the key.
   */
  public SortedMap<E, Integer> countMap() {
    final SortedMap<E, Integer> counts = new TreeMap<>();
    for (final Entry<E, Counter> entry : mMap.entrySet()) {
      counts.put(entry.getKey(), entry.getValue().count());
    }
    return counts;
  }

}
