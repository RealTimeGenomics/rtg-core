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

import java.util.Arrays;

import junit.framework.TestCase;

/**
 * Test Class
 */
public class QuickSortTest extends TestCase {


  public void testSort() {
    PortableRandom pr = new PortableRandom();
    final long seed = pr.getSeed();
    try {
      final int[] unsorted = new int[100000];
      for (int i = 0; i < unsorted.length; i++) {
        unsorted[i] = Math.abs(pr.nextInt());
      }
      final int[] unsorted2 = Arrays.copyOf(unsorted, unsorted.length);
      QuickSort.SortProxy proxy = new QuickSort.SortProxy() {

        @Override
        public int compare(long index1, long index2) {
          return Integer.valueOf(unsorted[(int) index1]).compareTo(unsorted[(int) index2]);
        }

        @Override
        public void swap(long index1, long index2) {
          int temp = unsorted[(int) index1];
          unsorted[(int) index1] = unsorted[(int) index2];
          unsorted[(int) index2] = temp;
        }

        @Override
        public long length() {
          return unsorted.length;
        }

        @Override
        public String toString() {
          return Arrays.toString(unsorted);

        }

      };
      QuickSort.sort(proxy);
      Arrays.sort(unsorted2);
      assertTrue(Arrays.equals(unsorted2, unsorted));
      assertTrue(QuickSort.isSorted(proxy));
    } catch (RuntimeException e) {
      throw new RuntimeException("Seed: " + seed + " has discovered an error", e);
    }
  }

}
