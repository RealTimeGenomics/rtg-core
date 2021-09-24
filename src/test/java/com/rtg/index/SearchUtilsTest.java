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
package com.rtg.index;

import java.util.Arrays;

import com.rtg.util.PortableRandom;
import com.rtg.util.array.CommonIndex;

import junit.framework.TestCase;

/**
 * Test class for search utils.
 */
public class SearchUtilsTest extends TestCase {

  public SearchUtilsTest(final String name) {
    super(name);
  }

  public void testBinarySearchLongStatic() {
    final CommonIndex a = dynamic(generateRandomArray(10000));
    int j = 0;
    for (int i = (int) a.get(0) - 1; i <= a.get(a.length() - 1) + 1; ++i) {
      if (j < a.length() && a.get(j) == i) {
        assertEquals(a.get(j), a.get(SearchUtils.binarySearch(a, 0, a.length() - 1, i)));
        while (j < a.length() && a.get(j) == i) {
          ++j;
        }
      } else {
        assertTrue(SearchUtils.binarySearch(a, 0, a.length() - 1, i) < 0);
      }
    }
  }

  public void testInterpolateSearchLongStatic() {
    final CommonIndex a = dynamic(generateRandomArray(10000));
    int j = 0;
    for (int i = (int) a.get(0) - 1; i <= a.get(a.length() - 1) + 1; ++i) {
      if (j < a.length() && a.get(j) == i) {
        assertEquals("i=" + i + " j=" + j, a.get(j), a.get(SearchUtils.interpolationSearch(a, 0, a.length() - 1, i)));
        while (j < a.length() && a.get(j) == i) {
          ++j;
        }
      } else {
        assertTrue(SearchUtils.interpolationSearch(a, 0, a.length() - 1, i) < 0);
      }
    }
  }

  public void testInterpolateSearch() {
    final int[] ra = generateRandomArray(100);
    //    for (int i = 0; i < ra.length; ++i) {
    //      if (i % 10 == 0) {
    //        System.err.println();
    //        System.err.print("[" + i + "] ");
    //      }
    //      System.err.print(" " + ra[i]);
    //    }
    //    System.err.println();
    final CommonIndex a = dynamic(ra);
    int j = 0;
    for (int i = (int) a.get(0) - 1; i <= a.get(a.length() - 1) + 1; ++i) {
      if (j < a.length() && a.get(j) == i) {
        assertEquals("i=" + i + " j=" + j, a.get(j), a.get(SearchUtils.interpolationSearch(a, 0, a.length() - 1, i)));
        while (j < a.length() && a.get(j) == i) {
          ++j;
        }
      } else {
        final long res = SearchUtils.interpolationSearch(a, 0, a.length() - 1, i);
        //System.err.println("res=" + res);
        assertTrue(res < 0);
      }
    }
  }

  public void testBracketSearch() {
    checkBracketSearch(1);
    checkBracketSearch(100);
    checkBracketSearch(10000);
  }

  private void checkBracketSearch(final int len) {
    final int[] ra = generateRandomArray(len);
    //    for (int i = 0; i < ra.length; ++i) {
    //      if (i % 10 == 0) {
    //        System.err.println();
    //        System.err.print("[" + i + "] ");
    //      }
    //      System.err.print(" " + ra[i]);
    //    }
    //    System.err.println();
    final CommonIndex a = dynamic(ra);
    for (long i = a.get(0); i < a.get(a.length() - 1); ++i) {
      final long index = SearchUtils.bracketSearch(a, 0, a.length() - 1, i);
      assertTrue(index >= 0);
      final long lowf = a.get(index);
      final long highf = a.get(index + 1);
      assertTrue(lowf <= i && i < highf);
    }
    assertEquals(-1, SearchUtils.bracketSearch(a, 0, a.length() - 1, -1));
    assertEquals(-1, SearchUtils.bracketSearch(a, 0, a.length() - 1, a.get(a.length() - 1)));
  }

  private static int[] generateRandomArray(final int len) {
    final int[] x = new int[len];
    final PortableRandom r = new PortableRandom(1);
    for (int i = 0; i < x.length; ++i) {
      if (i % 3 == 0) {
        x[i] = r.nextInt(20);
      } else {
        x[i] = r.nextInt(10000);
      }
    }
    Arrays.sort(x);
    return x;
  }

  private static CommonIndex dynamic(final int[] x) {
    final CommonIndex a = com.rtg.util.array.intindex.IntCreate.createIndex(x.length);
    for (int i = 0; i < x.length; ++i) {
      a.set(i, x[i]);
    }
    return a;
  }
}
