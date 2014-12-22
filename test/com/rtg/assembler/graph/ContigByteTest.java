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

package com.rtg.assembler.graph;

import junit.framework.TestCase;

/**
 */
public class ContigByteTest extends TestCase {
  public void test() {
    ContigByte contig = new ContigByte(new byte[] {1, 1, 2, 3, 3, 0, 4, 1});
    assertEquals(1, contig.nt(0));
    assertEquals(4, contig.nt(6));
    assertEquals(1, contig.nt(7));
    assertEquals(8, contig.length());
  }
}
