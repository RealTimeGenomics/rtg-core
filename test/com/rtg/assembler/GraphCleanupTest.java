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
    for (long i = 1; i <= graph.numberContigs(); i++) {
      assertEquals(isDeleted[(int) i - 1], graph.contigDeleted(i));
    }
  }
}
