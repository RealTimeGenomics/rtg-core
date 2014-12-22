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
      , Collections.<String, String>emptyMap()
      , Collections.<String, String>emptyMap()
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
