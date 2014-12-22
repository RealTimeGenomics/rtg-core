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
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * An unordered multi-set
 * @param <E> type of elements in the multi-set.
 */
public class MultiSet<E> extends IntegralAbstract {

  protected final Map<E, Counter> mMap;

  private long mTotalCount = 0;

  /**
   * Constructor
   */
  public MultiSet() {
    this(new HashMap<E, Counter>());
  }

  /**
   * Constructor with custom backing map
   * @param map the user-supplied map.
   */
  public MultiSet(Map<E, Counter> map) {
    mMap = map;
  }

  /**
   * Add a new element to the multi-set and increment its count.
   * @param element to be added.
   * @return count after add has been completed (&ge; 1).
   */
  public int add(E element) {
    final Counter c0 = mMap.get(element);
    final Counter c;
    if (c0 == null) {
      c = new Counter();
      mMap.put(element, c);
    } else {
      c = c0;
    }
    c.increment();
    mTotalCount++;
    return c.count();
  }

  /**
   * Add a new element to the multi-set and increment its count.
   * @param element to be added.
   * @param count to increment by.
   * @return count after add has been completed (&ge; 1).
   */
  public int add(final E element, final int count) {
    assert count >= 0;
    final Counter c0 = mMap.get(element);
    final Counter c;
    if (c0 == null) {
      if (count == 0) {
        return 0;
      }
      c = new Counter();
      mMap.put(element, c);
    } else {
      c = c0;
    }
    c.increment(count);
    mTotalCount += count;
    return c.count();
  }

  /**
   * Get the count for the element.
   * @param element to be fetched.
   * @return the count for the element (0 if not in the multi-set otherwise &ge; 1).
   */
  public int get(E element) {
    final Counter c0 = mMap.get(element);
    if (c0 == null) {
      return 0;
    } else {
      return c0.count();
    }
  }

  /**
   * Get the set of keys in this map
   * @return the set of keys
   */
  public Set<E> keySet() {
    return mMap.keySet();
  }

  /**
   * @return total count
   */
  public long totalCount() {
    return mTotalCount;
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append("[ ");
    int count = 0;
    int lines = 1;
    for (final Entry<E, Counter> entry : mMap.entrySet()) {
      if (count > 0) {
        if (count % 10 == 0) {
          sb.append(StringUtils.LS);
          lines++;
        }
        sb.append(", ");
      }
      sb.append(entry.getKey().toString()).append("->").append(entry.getValue().count());
      count++;
    }
    sb.append(lines > 1 ? StringUtils.LS : "").append("]");
    return sb.toString();
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    long tot = 0;
    for (final E key : keySet()) {
      tot += get(key);
    }
    Exam.assertEquals(tot, mTotalCount);
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mTotalCount >= mMap.size());
    return true;
  }

}
