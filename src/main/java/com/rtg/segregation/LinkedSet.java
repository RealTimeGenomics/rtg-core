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

package com.rtg.segregation;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NavigableMap;
import java.util.TreeMap;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Maintains a set of elements and an ability to iterate over all elements between two specified elements.
 */
class LinkedSet<E> extends IntegralAbstract {

  private final NavigableMap<Integer, E> mElements = new TreeMap<>();
  private final Map<E, Integer> mIndexes = new HashMap<>();

  boolean add(E element) {
    if (mIndexes.containsKey(element)) {
      return false;
    }
    final int size = mElements.size();
    mElements.put(size, element);
    mIndexes.put(element, size);
    return true;
  }

  Iterator<E> iterator(E start, E end) {
    final Integer startIndex = mIndexes.get(start);
    final Integer endIndex = mIndexes.get(end);
    return mElements.subMap(startIndex, false, endIndex, true).values().iterator();
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (final Entry<E, Integer> el : mIndexes.entrySet()) {
      final E key = el.getKey();
      final Integer value = el.getValue();
      Exam.assertTrue(mElements.get(value) == key);
    }
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mElements);
    Exam.assertNotNull(mIndexes);
    Exam.assertEquals(mElements.size(), mIndexes.size());
    return true;
  }

}
