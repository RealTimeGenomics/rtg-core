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
package com.rtg.ngs.blocking;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 * Tests corresponding class
 */
public class MapQScoringReadBlockerTest extends TestCase {

  protected static final int READ_ID0 = 0;
  protected static final int READ_ID1 = (1 << 16) - 1; // a bigger readid.

  protected String expectedDiagnosticString() {
    return "blocked pairings";
  }

  protected String expectedToString() {
    return "ScoringReadBlocker";
  }

  protected MapQScoringReadBlocker getScoringReadBlocker(final int count, final int threshold) {
    return new MapQScoringReadBlocker(count, threshold, expectedDiagnosticString());
  }

  public void testNorm() {
    final MapQScoringReadBlocker srb = new MapQScoringReadBlocker(1, 20, "my title");
    assertEquals(23, srb.getMapQ(0));
    assertEquals(0, srb.getCount1(0));
    srb.increment(0, 0); // First good hit
    assertEquals(1, srb.getCount1(0));
    int lmq = srb.getMapQ(0);
    assertEquals(37, lmq);

    assertEquals(0, srb.getCount2(0));
    for (int i = 0; i < 40; i++) {  // Fire in some second best hits
      srb.increment(0, 5);
      final int mq = srb.getMapQ(0);
      //System.out.println("i=" + i + " mq " + mq);
      assertTrue(mq <= lmq);        // Should be decreasing with additional second best hits
      lmq = mq;
    }
    assertEquals(40, srb.getCount2(0));
    for (int i = 0; i < 19; i++) { // Some more at the best score
      srb.increment(0, 0);
    }
    assertEquals(0, srb.getMapQ(0)); // Multiple at best core gets mapq=0
  //srb.dumpRead(0);
    assertFalse(srb.isBlocked1(0, 0));
    assertFalse(srb.isBlocked2(0, 5));
    srb.increment(0, 0);
  //srb.dumpRead(0);
    assertTrue(srb.isBlocked1(0, 0));
    assertFalse(srb.isBlocked2(0, 5));
    for (int i = 0; i < 255; i++) {
      srb.increment(0, 5);
    }
    //srb.dumpRead(0);
    assertTrue(srb.isBlocked2(0, 5));
  }

  public void testEdge() {
    final MapQScoringReadBlocker srb = new MapQScoringReadBlocker(1, 255);
    for (int i = 0; i < 255; i++) {
      srb.increment(0, 0);
    }
    //srb.dumpRead(0);
    assertFalse(srb.isBlocked1(0, 0));
    assertFalse(srb.isBlocked2(0, 5));
    srb.increment(0, 0);
    //srb.dumpRead(0);
    assertTrue(srb.isBlocked1(0, 0));
    assertFalse(srb.isBlocked2(0, 5));
  }

  public void testConstant() {
    for (int k = 1; k < 256; k <<= 1) {
      final int count = READ_ID1 + 1;
      final MapQScoringReadBlocker b = getScoringReadBlocker(count, k);
      for (int j = 0; j <= k; j++) {
        assertFalse(b.isBlocked1(READ_ID1, (byte) 1));
        assertFalse(b.isBlocked1(READ_ID0, (byte) 1));
        b.increment(READ_ID1, (byte) 1);
      }
      assertTrue(b.isBlocked1(READ_ID1, (byte) 1));
      assertFalse(b.isBlocked1(READ_ID0, (byte) 1));
    }
  }

  public void testIncreasing() {
    for (int k = 1; k < 256; k <<= 1) {
      final MapQScoringReadBlocker b = getScoringReadBlocker(2, k);
      assertEquals(expectedToString(), b.toString());
      assertFalse(b.isBlocked1(1, 0));
      assertFalse(b.isBlocked1(0, 0));
      b.increment(1, (byte) 0);
      for (int j = 1; j <= k; j++) {
        final int score = k * 2;
        assertTrue(b.isBlocked1(1, score));
        assertTrue(b.isBlocked1(0, score) == (score > 255));
        b.increment(1, 0);
      }
      assertTrue(b.isBlocked1(1, 0));
      assertFalse(b.isBlocked1(0, 0));
      final MemoryPrintStream ps = new MemoryPrintStream();
      Diagnostic.setLogStream(ps.printStream());
      assertEquals(0, b.getCount1(0));
      assertEquals(255, b.getScore1(0));
      b.close();
      assertTrue(ps.toString(), ps.toString().contains(expectedDiagnosticString()));
      Diagnostic.setLogStream();
      ps.close();
    }
  }

  public void testDecreasing() {
    final MapQScoringReadBlocker b = getScoringReadBlocker(1, 10);
    for (int score = 512; score > 0; score >>= 1) {
      assertFalse(b.isBlocked1(0, 0));
      for (int j = 0; j < 12; j++) {
        //System.err.println("score=" + k + " j=" + j);
        assertTrue(b.isBlocked1(0, score) == ((j > 10) || (score > 255)));
        assertTrue(b.isBlocked2(0, score) == (score > 255));
        b.increment(0, score);
      }
      assertFalse(b.isBlocked1(0, 0));
    }
    assertEquals(1, b.getScore1(0));
    assertEquals(12, b.getCount1(0));
    assertEquals(2, b.getScore2(0));
    assertEquals(12, b.getCount2(0));
  }

  public void testCons() {
    try {
      getScoringReadBlocker(0, 0);
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals(e.getMessage(), "threshold must be > 0, not 0"); // ok
    }
    try {
      getScoringReadBlocker(0, 65536);
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals(e.getMessage(), "threshold must be > 0, not 65536"); // ok
    }
  }

}
