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
  List<Integer> ref(int ... refInts) {
    final List<Integer> res = new ArrayList<>();
    for (int i : refInts) {
      res.add(i);
    }
    return res;
  }
}
