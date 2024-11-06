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
package com.rtg.util;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * A map from a key to a set of values.
 *
 * @param <K> type of keys
 * @param <V> type of values
 */
public class MapSet<K, V> extends IntegralAbstract {

  private final Map<K, Set<V>> mMap = new HashMap<>();

  /**
   * Add a new key value pair.
   *
   * @param key index to be added.
   * @param value associated value to be added.
   */
  public void put(final K key, final V value) {
    final Set<V> set = mMap.get(key);
    final Set<V> nset;
    if (set == null) {
      nset = new HashSet<>();
      mMap.put(key, nset);
    } else {
      nset = set;
    }
    nset.add(value);
  }

  /**
   * Remove a key value pair.
   *
   * @param key used for search and then removal.
   * @return the previous set if any (null otherwise)
   */
  public Set<V> remove(final K key) {
    return mMap.remove(key);
  }

  /**
   * Get the key set.
   * @return the keyset
   */
  public Set<K> keySet() {
    return mMap.keySet();
  }

  /**
   * Check if key value pair is in the map.
   * @param key to be used for search.
   * @param value expected to be found.
   * @return true iff the pair is in the map.
   */
  public boolean contains(final K key, final V value) {
    final Set<V> set = mMap.get(key);
    return set != null && set.contains(value);
  }

  /**
   * Get the set of values associated with the key.
   * @param key to be used for search.
   * @return the set of values associated with the key - null if there are no
   * values stored.
   */
  public Set<V> get(final K key) {
    return mMap.get(key);
  }

  /**
   * Get a set of values in the map.
   * @return the set of values in the map.
   */
  public Set<Map.Entry<K, Set<V>>> entrySet() {
    return mMap.entrySet();
  }

  /**
   * Get the number of different keys.
   * @return the number of different keys.
   */
  public int numberOfKeys() {
    return mMap.keySet().size();
  }

  @Override
  public int hashCode() {
    return mMap.hashCode() + 31;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    if (getClass() != obj.getClass()) {
      return false;
    }
    final MapSet<?, ?> other = (MapSet<?, ?>) obj;
    return mMap.equals(other.mMap);
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mMap);
    return true;
  }

  @Override
  public String toString() {
    final Map<K, Set<V>> sorted = new TreeMap<>(mMap);
    return sorted.toString();
  }
}

