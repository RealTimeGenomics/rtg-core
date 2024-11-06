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
