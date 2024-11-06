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

import java.io.IOException;
import java.util.Collections;

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
    final GraphKmerAttribute a = GraphMapCliTest.makeGraph(4, new String[] {"AAAACCCC", "TTTTTGG", "CCCCC", "ATATA"}, new long[][] {{1, 2}, {2, 3}, {3, 4}});
    final GraphKmerAttribute b = GraphMapCliTest.makeGraph(4, new String[] {"TTTTTGG", "AAAACCCC", "GGGGG", "TATAT"}, new long[][] {{2, 1}, {1, -3}, {-3, -4}});
    final Graph sortedA = GraphSorter.sortedGraph(a);
    final Graph sortedB = GraphSorter.sortedGraph(b);
    final StoreDirString aString = new StoreDirString();
    final StoreDirString bString = new StoreDirString();
    GraphWriter.write(sortedA, aString, "foo", Collections.emptySet());
    GraphWriter.write(sortedB, bString, "foo", Collections.emptySet());
    final String aFiltered = StringUtils.grepMinusV(aString.toString(), "header.tsv|guid|date");
    final String bFiltered = StringUtils.grepMinusV(bString.toString(), "header.tsv|date|guid");
    assertEquals(aFiltered, bFiltered);
  }
}
