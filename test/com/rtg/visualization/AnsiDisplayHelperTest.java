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

package com.rtg.visualization;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

/**
 */
public class AnsiDisplayHelperTest {

  @Test
  public void testAnsiColors() {
    assertEquals((char) 27 + "[48;5;17m", AnsiDisplayHelper.ansiBackground(DisplayHelper.BLUE));
    assertEquals((char) 27 + "[31m", AnsiDisplayHelper.ansiForeground(DisplayHelper.RED));
  }

  @Rule
  public ExpectedException expectedException = ExpectedException.none();
  @Test
  public void testColorExceptionRed() {
    expectedException.expect(IllegalArgumentException.class);
    AnsiDisplayHelper.extendedColor(6, 5, 5);
  }
  @Test
  public void testColorExceptionGreen() {
    expectedException.expect(IllegalArgumentException.class);
    AnsiDisplayHelper.extendedColor(5, 6, 5);
  }
  @Test
  public void testColorExceptionBlue() {
    expectedException.expect(IllegalArgumentException.class);
    AnsiDisplayHelper.extendedColor(5, 5, 6);
  }
  @Test
  public void testExtendedColorWhite() {
    assertEquals(231, AnsiDisplayHelper.extendedColor(5, 5, 5));
  }

  @Test
  public void testMarkupStart() {
    final AnsiDisplayHelper helper = new AnsiDisplayHelper();
    assertTrue(helper.isMarkupStart((char) 0x1b));
    assertFalse(helper.isMarkupStart('f'));
  }

  @Test
  public void testMarkupEnd() {
    final AnsiDisplayHelper helper = new AnsiDisplayHelper();
    assertTrue(helper.isMarkupEnd('m'));
    assertFalse(helper.isMarkupEnd((char) 0x1b));
    assertFalse(helper.isMarkupEnd((char) 0x00));
  }
}
