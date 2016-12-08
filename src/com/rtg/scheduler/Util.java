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

package com.rtg.scheduler;

import java.util.Collection;

import com.rtg.util.integrity.Exam;

/**
 */
public final class Util {
  private Util() { }

  /**Check that x is ordered with respect to everything in the set.
   * @param x item to be checked.
   * @param set set to be checked.
   * @param order ordering to be used.
   * @param <X> type of items being compared.
   * @return true iff x is ordered consistently with everything in set.
   */
  public static <X> boolean checkOrder(final Comparable<X> x, Collection<X> set, final int order) {
    //System.err.println("x=" + x + " order=" + order + " set=" + set);
    for (final X y : set) {
      if (y != null) {
        final int c = x.compareTo(y);
        Exam.assertTrue(/*"  y=" + y + " compare=" + x.compareTo(y), */order < 0 ? c > 0 : c < 0);
      }
    }
    return true;
  }

  /**
   * Count the number of non-null objects in the collection.
   * @param collection to be counted.
   * @return the number of non-null objects. in collection.
   */
  public static int nonNullSize(final Collection<?> collection) {
    int count = 0;
    for (final Object obj : collection) {
      if (obj != null) {
        ++count;
      }
    }
    return count;
  }
}
