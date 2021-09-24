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

import java.util.TreeMap;

import com.rtg.assembler.graph.Graph;

import junit.framework.TestCase;

/**
 */
public class ContigPositionTest extends TestCase {
  public void testGraphPosition() {
    final Graph g = GraphMapCliTest.makeGraph(0, new String[]{"AAAAA", "GGGGG", "CCC"}, new long[][]{{1, 2}});
    final TreeMap<Long, Long> positionDecoder = ContigPosition.buildDecoder(g);
    assertEquals(0, new ContigPosition(1, 0, g).encode());
    assertEquals(1, new ContigPosition(1, 1, g).encode());
    assertEquals(5, new ContigPosition(2, 0, g).encode());
    assertEquals(9, new ContigPosition(2, 4, g).encode());
    assertEquals(10, new ContigPosition(3, 0, g).encode());
    ContigPosition contigPosition;
    contigPosition = new ContigPosition(0, g, positionDecoder);
    assertEquals(1, contigPosition.mContigId);
    assertEquals(0, contigPosition.mPosition);

    contigPosition = new ContigPosition(9, g, positionDecoder);
    assertEquals(2, contigPosition.mContigId);
    assertEquals(4, contigPosition.mPosition);

    contigPosition = new ContigPosition(10, g, positionDecoder);
    assertEquals(3, contigPosition.mContigId);
    assertEquals(0, contigPosition.mPosition);

  }
}
