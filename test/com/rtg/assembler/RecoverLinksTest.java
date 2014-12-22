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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;

import junit.framework.TestCase;

/**
 */
public class RecoverLinksTest extends TestCase {
  public void testLink() {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(5, new String[]{"AAAACCACCACAGTATTAGATG", "AGGGGTGGGTGAGACC", "ATTTTGTTGTG", "TTGGGAGATTGA"}, new long[][]{});
    ConstraintCache cache = new ConstraintCache();
    cache.addConstraint(1, -2, 3, 4, 20);
    cache.addConstraint(-2, 1, 7, 12, 20);

    cache.addConstraint(2, -3, 7, 9, 20);
    cache.addConstraint(2, -4, 7, 9, 20);
    RecoverLinks.recover(Arrays.asList(cache), graph);
    final List<Long> path = new ArrayList<>();
    for (int i = 0; i < graph.pathLength(1); i++) {
      path.add(graph.pathContig(1, i));
    }
    assertEquals(Arrays.asList(-2L, 5L, -1L), path);
    assertEquals("AGATGNNNNNNNAGGGG", ContigString.contigSequenceString(graph.contig(-5)));
    assertEquals(5, graph.numberContigs());
  }
}
