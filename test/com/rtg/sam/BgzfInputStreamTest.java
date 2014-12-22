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
package com.rtg.sam;

import com.rtg.util.Resources;

import junit.framework.TestCase;

/**
 * Tests corresponding class
 */
public class BgzfInputStreamTest extends TestCase {

  public void testInput() throws Exception {

    BgzfInputStream bgzfi = new BgzfInputStream(Resources.getResourceAsStream("com/rtg/sam/resources/bam.bam"));
    try {

      assertEquals(0L, bgzfi.blockStart());

      final byte[] b = new byte[197];
      final int i = bgzfi.read(b);
      assertEquals(197, i);
      assertEquals(0L, bgzfi.blockStart());

    } finally {
      bgzfi.close();
    }
    bgzfi = new BgzfInputStream(Resources.getResourceAsStream("com/rtg/sam/resources/bam.bam"));
    try {
      final byte[] b = new byte[199];

      try {
        final int len = bgzfi.read(b, 10, 199);
        fail();
        assertEquals("wasteful line of code this", 0, len);
      } catch (ArrayIndexOutOfBoundsException aioobe) {
      }

      final int i = bgzfi.read(b);
      assertEquals(198, i);
      assertEquals(134L, bgzfi.blockStart());

    } finally {
      bgzfi.close();
    }
  }
}
