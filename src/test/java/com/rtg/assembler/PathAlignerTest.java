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

import java.util.Collections;

import com.rtg.assembler.graph.MutableGraph;
import com.rtg.mode.DnaUtils;

import junit.framework.TestCase;


/**
 */
public class PathAlignerTest extends TestCase {

  private MutableGraph getGraph() {
    return GraphMapCliTest.makeGraph(1, new String[] {"AAAAAC", "CTAGTCA"}
      , new long[][] {{1, 2, 1, 2, 1, 2, 1, 2}}
      , Collections.emptyMap()
      , Collections.emptyMap()
      );
  }

  public void testExactMatch() {
    final PathAligner aligner = new PathAligner(getGraph());
    //                   AAAAACTAGTCAAAAACTAGTCAAAAACTAGTCAAAAACTAGTCA
    final String read = "AAAAACTAGTC";
    final byte[] frag = DnaUtils.encodeString(read);
    assertEquals(0, aligner.align(frag, 0, 1, 0));
  }

  public void testLongerExactMatch() {
    final PathAligner aligner = new PathAligner(getGraph());
    //                   AAAAACTAGTCAAAAACTAGTCAAAAACTAGTCAAAAACTAGTCA
    final String read = "AAAAACTAGTCAAAAACTAGT";
    final byte[] frag = DnaUtils.encodeString(read);
    assertEquals(0, aligner.align(frag, 0, 1, 0));
  }

  public void testOneInsertion() {
    final PathAligner aligner = new PathAligner(getGraph());
    //                   AAAAACTA-GTCAAAAACTAGTCAAAAACTAGTCAAAAACTAGTCA
    final String read = "AAAAACTATGTCAAAAACTAGT";
    final byte[] frag = DnaUtils.encodeString(read);
    assertEquals(2, aligner.align(frag, 0, 1, 0));
  }

  public void testOneDeletion() {
    final PathAligner aligner = new PathAligner(getGraph());
    //                   AAAAACTAGTCAAAAACTAGTCAAAAACTAGTCAAAAACTAGTCA
    final String read = "AAAAACTA-TCAAAAACTAGT";
    final byte[] frag = DnaUtils.encodeStringWithHyphen(read);
    assertEquals(2, aligner.align(frag, 0, 1, 0));
  }

  public void testMixed() {
    final PathAligner aligner = new PathAligner(getGraph());
    //                   AAAAACTAGTCAAAAACTAGTCAAAAACTAG-TCAAA-AACTAGTCA
    final String read = "AAAAACTA-TCAAAAACTAGTCAAA--CTAGCTCAAATAACTAG";
    final byte[] frag = DnaUtils.encodeStringWithHyphen(read);
    assertEquals(9, aligner.align(frag, 0, 1, 0));
  }

}
