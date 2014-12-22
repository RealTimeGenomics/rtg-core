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

import java.util.TreeMap;

/**
 * Map that will automatically add a default value as produced by the make method if the requested key isn't currently
 * there.
 * Should <code>getOrAdd</code> override regular get method?
 */
public abstract class AutoAddMap<K, V> extends TreeMap<K, V> {
  /**
   * Retrieve the value stored at <code>name</code> or a new instance if it doesn't exist
   * @param name key of value to retrieve
   * @return the stored value or a new one
   */
  public V getOrAdd(K name) {
    if (containsKey(name)) {
      return get(name);
    }
    final V res = make();
    put(name, res);
    return res;
  }

  /**
   * @return a newly constructed instance of the stored type
   */
  public abstract V make();
}
