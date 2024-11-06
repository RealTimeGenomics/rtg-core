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
