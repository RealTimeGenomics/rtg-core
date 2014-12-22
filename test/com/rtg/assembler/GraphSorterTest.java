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

import java.io.IOException;
import java.util.Collections;
import java.util.UUID;

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.io.GraphWriter;
import com.rtg.util.StringUtils;
import com.rtg.util.store.StoreDirString;

import junit.framework.TestCase;

/**
 */
public class GraphSorterTest extends TestCase {
  public void testGraphSorter() throws IOException {
    GraphKmerAttribute a = GraphMapCliTest.makeGraph(4, new String[] {"AAAACCCC", "TTTTTGG", "CCCCC", "ATATA"}, new long[][] {{1, 2}, {2, 3}, {3, 4}});
    GraphKmerAttribute b = GraphMapCliTest.makeGraph(4, new String[] {"TTTTTGG", "AAAACCCC", "GGGGG", "TATAT"}, new long[][] {{2, 1}, {1, -3}, {-3, -4}});
    Graph sortedA = GraphSorter.sortedGraph(a);
    Graph sortedB = GraphSorter.sortedGraph(b);
    StoreDirString aString = new StoreDirString();
    StoreDirString bString = new StoreDirString();
    GraphWriter.write(sortedA, aString, "foo", Collections.<UUID>emptySet());
    GraphWriter.write(sortedB, bString, "foo", Collections.<UUID>emptySet());
    String aFiltered = StringUtils.grepMinusV(aString.toString(), "header.tsv|guid|date");
    String bFiltered = StringUtils.grepMinusV(bString.toString(), "header.tsv|date|guid");
    assertEquals(aFiltered, bFiltered);

  }
}
