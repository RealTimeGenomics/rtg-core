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
