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
public class ProteinFrameTest extends TestCase {

  /**
   * Test method for {@link com.rtg.mode.ProteinFrame}.
   */
  public final void test() {
    TestUtils.testPseudoEnum(ProteinFrame.class, "[PROTEIN]");
    try {
      ProteinFrame.PROTEIN.getReverse();
      fail();
    } catch (UnsupportedOperationException e) {
      assertEquals("Not supported", e.getMessage());
    }
  }

  /**
   * Test method for {@link com.rtg.mode.ProteinFrame#display()}.
   */
  public final void testCode() {
    for (final ProteinFrame bi : ProteinFrame.values()) {
      assertEquals(bi, ProteinFrame.frameFromCode(bi.ordinal()));
    }
  }

  public final void testMode() {
    ProteinFrame.values();
  }

  /**
   * Test method for {@link com.rtg.mode.ProteinFrame#phase()}.
   */
  public final void testPhase() {
    for (final ProteinFrame pr : ProteinFrame.values()) {
      assertEquals(0, pr.phase());
    }
  }

  /**
   * Test method for {@link com.rtg.mode.ProteinFrame#display()}.
   */
  public final void testDisplay() {
    assertEquals("", ProteinFrame.PROTEIN.display());
  }

  /**
   * Test method for {@link com.rtg.mode.ProteinFrame#isForward()}.
   */
  public final void testIsForward() {
    assertTrue(ProteinFrame.PROTEIN.isForward());
  }

  /**
   * Test method for {@link com.rtg.mode.ProteinFrame#frameFromCode(int)}.
   */
  public final void testFrameFromCodeBad() {
    try {
      ProteinFrame.frameFromCode(-1);
      fail("IllegalArgumentException expected");
    } catch (final IllegalArgumentException e) {
      assertEquals("-1", e.getMessage());
    }
    try {
      ProteinFrame.frameFromCode(1);
      fail("IllegalArgumentException expected");
    } catch (final IllegalArgumentException e) {
      assertEquals("1", e.getMessage());
    }
  }

  /**
   * Test method for {@link com.rtg.mode.ProteinFrame#code(byte[], int, int)}.
   */
  public final void testCode1() {
    final byte[] codes = {0, 1, 2, 3};
    final Frame f = ProteinFrame.PROTEIN;
    assertEquals(0, f.code(codes, 3, 0));
    assertEquals(1, f.code(codes, 3, 1));
    assertEquals(2, f.code(codes, 3, 2));
    assertEquals(0, f.code(codes, 3, 3));
    try {
      f.code(codes, 3, -1);
      fail("RuntimeException expected");
    } catch (final RuntimeException e) {
      //expected
    }
  }
}

