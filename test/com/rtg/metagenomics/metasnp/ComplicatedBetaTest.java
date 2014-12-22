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
package com.rtg.metagenomics.metasnp;

import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;

/**
 */
public class ComplicatedBetaTest extends TestCase {
  public void testInitZero() {
    final List<int[]> assignments = new ArrayList<>();
    assignments.add(new int[] {3, 3, 3});
    final ComplicatedBeta beta = new ComplicatedBeta(ref(3), assignments, 3);
    assertEquals(2, beta.getCount(2, new int[] {2, 2, 2}));
  }
  public void testInitSecondRef() {
    final List<int[]> assignments = new ArrayList<>();
    assignments.add(new int[] {3, 3, 3});
    final ComplicatedBeta beta = new ComplicatedBeta(ref(3), assignments, 3);
    assertEquals(2, beta.getCount(2, new int[] {2, 2, 2}));
  }
  public void testInitMax() {
    final List<int[]> assignments = new ArrayList<>();
    assignments.add(new int[] {0, 1, 2});
    assignments.add(new int[] {3, 2, 0});
    final ComplicatedBeta beta = new ComplicatedBeta(ref(3, 1), assignments, 3);
    assertEquals(3, beta.getCount(2, new int[] {1, 0, 3}));
  }
  public void testManyStrains() {
    final List<int[]> assignments = new ArrayList<>();
    assignments.add(new int[] {0, 1, 2, 3, 0, 1});
    assignments.add(new int[] {1, 0, 3, 2, 1, 0});
    final ComplicatedBeta beta = new ComplicatedBeta(ref(3, 2), assignments, 3);
    assertEquals(3, beta.getCount(1, new int[]{3, 2, 0, 1, 3, 2}));

  }

  public void testMulti() {
    final List<int[]> assignments = new ArrayList<>();
    assignments.add(new int[] {3, 0, 3});
    assignments.add(new int[] {0, 1, 1});
    assignments.add(new int[] {0, 1, 3});
    assignments.add(new int[] {2, 2, 2});
    final ComplicatedBeta beta = new ComplicatedBeta(ref(3, 0, 2, 1), assignments, 3);
    assertEquals(2, beta.getCount(2, new int[] {1, 1, 1}));
    assertEquals(1, beta.getCount(2, new int[]{1, 1, 2}));
  }
  List<Byte> ref(int ... refInts) {
    final List<Byte> res = new ArrayList<>();
    for (int i : refInts) {
      res.add((byte) i);
    }
    return res;
  }
}
