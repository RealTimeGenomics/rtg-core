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

package com.rtg.assembler;

import com.rtg.assembler.graph.implementation.GraphKmerAttribute;

import junit.framework.TestCase;

/**
 */
public class GraphCleanupTest extends TestCase {
  public void testClean() {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(4, new String[]{"AAAAA", "TTATTATATA", "ACGTG", "ATCACGTGTAGATAG", "AAAATAAATAACA", "ACGTA", "ACGTG", "ACCGT", "AGGTG", "AGGTG", "ACGAG"}, new long[][]{{2, 3}, {3, 4}, {6, 6}, {8, 7}, {7, 9}, {4, 10}});
    final int deleted = GraphCleanup.clean(5, graph);
    assertEquals(7, deleted);
    final boolean[]  isDeleted = {true, false, false, false, false, true, true, true, true, true, true};
    assertEquals(isDeleted.length, graph.numberContigs());
    for (long i = 1; i <= graph.numberContigs(); ++i) {
      assertEquals(isDeleted[(int) i - 1], graph.contigDeleted(i));
    }
  }
}
