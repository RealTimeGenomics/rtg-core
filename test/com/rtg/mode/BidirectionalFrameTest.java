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
package com.rtg.mode;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class BidirectionalFrameTest extends TestCase {

  /**
   * Test method for {@link com.rtg.mode.BidirectionalFrame}.
   */
  public final void test() {
    TestUtils.testPseudoEnum(BidirectionalFrame.class, "[FORWARD, REVERSE]");
    assertEquals(BidirectionalFrame.FORWARD, BidirectionalFrame.REVERSE.getReverse());
    assertEquals(BidirectionalFrame.REVERSE, BidirectionalFrame.FORWARD.getReverse());
  }

  /**
   * Test method for {@link com.rtg.mode.BidirectionalFrame}.
   */
  public final void testValues() {
    for (final BidirectionalFrame bi : BidirectionalFrame.values()) {
      assertEquals(bi, BidirectionalFrame.frameFromOrdinal(bi.ordinal()));
    }
  }

  /**
   * Test method for {@link com.rtg.mode.BidirectionalFrame#display()}.
   */
  public final void testDisplay() {
    assertEquals("F", BidirectionalFrame.FORWARD.display());
    assertEquals("R", BidirectionalFrame.REVERSE.display());
  }

  /**
   * Test method for {@link com.rtg.mode.BidirectionalFrame#display()}.
   */
  public final void testIsForward() {
    assertTrue(BidirectionalFrame.FORWARD.isForward());
    assertFalse(BidirectionalFrame.REVERSE.isForward());
  }

  /**
   * Test method for {@link com.rtg.mode.BidirectionalFrame#phase()}.
   */
  public final void testPhase() {
    for (final BidirectionalFrame bi : BidirectionalFrame.values()) {
      assertEquals(0, bi.phase());
    }
  }

  /**
   * Test method for {@link com.rtg.mode.TranslatedFrame#frameFromCode(int)}.
   */
  public final void testFrameFromCodeBad() {
    try {
      BidirectionalFrame.frameFromOrdinal(-1);
      fail("Exception expected");
    } catch (final RuntimeException e) {
      //expected
    }
    try {
      BidirectionalFrame.frameFromOrdinal(2);
      fail("RuntimeException expected");
    } catch (final RuntimeException e) {
      //expected
    }
  }

  /**
   * Test method for {@link com.rtg.mode.TranslatedFrame#code(byte[], int, int)}.
   */
  public final void testCodeF() {
    final byte[] codes = {0, 1, 2, 3};
    final Frame f = BidirectionalFrame.FORWARD;
    assertEquals(0, f.code(codes, 3, 0));
    assertEquals(1, f.code(codes, 3, 1));
    assertEquals(2, f.code(codes, 3, 2));
    try {
      f.code(codes, 3, 3);
      fail("RuntimeException expected");
    } catch (final RuntimeException e) {
      assertEquals("length=3 index=3", e.getMessage());
    }
    try {
      f.code(codes, 3, -1);
      fail("RuntimeException expected");
    } catch (final RuntimeException e) {
      //expected
    }
  }

  public final void testCodeR() {
    final byte[] codes = {0, 1, 2, 3, 4};
    final Frame f = BidirectionalFrame.REVERSE;
    assertEquals(1, f.code(codes, 5, 0));
    assertEquals(2, f.code(codes, 4, 0));
    assertEquals(3, f.code(codes, 3, 0));
    assertEquals(4, f.code(codes, 3, 1));
    assertEquals(0, f.code(codes, 3, 2));
    try {
      f.code(codes, 3, 3);
      fail("RuntimeException expected");
    } catch (final RuntimeException e) {
      //expected
    }
    try {
      f.code(codes, 3, -1);
      fail("RuntimeException expected");
    } catch (final RuntimeException e) {
      assertEquals("length=3 index=-1", e.getMessage());
    }
  }

  public final void testComplement() {
    for (final DNA dna : DNA.values()) {
      assertEquals(dna.complement(), DNA.values()[BidirectionalFrame.complement((byte) dna.ordinal())]);
    }
  }
}

