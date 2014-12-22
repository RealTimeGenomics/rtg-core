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
public class UnidirectionalFrameTest extends TestCase {

  /**
   * Test method for {@link com.rtg.mode.UnidirectionalFrame}.
   */
  public final void test() {
    TestUtils.testPseudoEnum(UnidirectionalFrame.class, "[FORWARD]");
    try {
      UnidirectionalFrame.FORWARD.getReverse();
      fail();
    } catch (UnsupportedOperationException e) {
      assertEquals("Not supported", e.getMessage());
    }
  }

  /**
   * Test method for {@link com.rtg.mode.UnidirectionalFrame#display()}.
   */
  public final void testCode() {
    for (final UnidirectionalFrame bi : UnidirectionalFrame.values()) {
      assertEquals(bi, UnidirectionalFrame.frameFromCode(bi.ordinal()));
    }
  }

  /**
   * Test method for {@link com.rtg.mode.UnidirectionalFrame#phase()}.
   */
  public final void testPhase() {
    for (final UnidirectionalFrame un : UnidirectionalFrame.values()) {
      assertEquals(0, un.phase());
    }
  }


  /**
   * Test method for {@link com.rtg.mode.UnidirectionalFrame#display()}.
   */
  public final void testDisplay() {
    assertEquals("", UnidirectionalFrame.FORWARD.display());
  }

  /**
   * Test method for {@link com.rtg.mode.UnidirectionalFrame#isForward()}.
   */
  public final void testIsForward() {
    assertTrue(UnidirectionalFrame.FORWARD.isForward());
  }

  /**
   * Test method for {@link com.rtg.mode.UnidirectionalFrame#frameFromCode(int)}.
   */
  public final void testFrameFromCodeBad() {
    final String vmName = System.getProperty("java.vm.name", "");
    try {
      UnidirectionalFrame.frameFromCode(-1);
      fail("RuntimeException expected");
    } catch (final RuntimeException e) {
      final String expected = vmName.startsWith("IBM J9") ? "Array index out of range: -1" : "-1";
      assertEquals(expected, e.getMessage());
    }
    try {
      UnidirectionalFrame.frameFromCode(1);
      fail("RuntimeException expected");
    } catch (final RuntimeException e) {
      final String expected = vmName.startsWith("IBM J9") ? "Array index out of range: 1" : "1";
      assertEquals(expected, e.getMessage());
    }
  }

  /**
   * Test method for {@link com.rtg.mode.ProteinFrame#code(byte[], int, int)}.
   */
  public final void testCode1() {
    final byte[] codes = {0, 1, 2, 3};
    final Frame f = UnidirectionalFrame.FORWARD;
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
}

