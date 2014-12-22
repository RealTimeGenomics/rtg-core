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
    Set<K> ret;
    //
    ret = mMap.keySet();
    return ret;
  }

  /**
   * Check if key value pair is in the map.
   * @param key to be used for search.
   * @param value expected to be found.
   * @return true iff the pair is in the map.
   */
  public boolean contains(final K key, final V value) {
    final Set<V> set = mMap.get(key);
    if (set == null) {
      return false;
    }
    return set.contains(value);
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

