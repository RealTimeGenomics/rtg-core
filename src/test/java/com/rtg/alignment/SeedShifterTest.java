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

package com.rtg.alignment;

import com.rtg.util.Utils;

import junit.framework.TestCase;

/**
 */
public class SeedShifterTest extends TestCase {

  public void test() {
    final SeedShifter ss = new SeedShifter(new Seed(3), new byte[] {1, 2, 3, 0, 1, 2, 3}, 7, -2);
    ss.integrity();
    assertEquals("position=-2 value=11111100:00000000:00000000:00111111", ss.toString());
    assertFalse(ss.isValid());

    assertEquals("11111100:00000000:00000000:00111111", Utils.toBitsSep(ss.next()));
    assertEquals("position=-1 value=11111100:00000000:00000000:00111111", ss.toString());
    assertFalse(ss.isValid());

    assertEquals("11111100:00000000:00000000:00111111", Utils.toBitsSep(ss.next()));
    assertEquals("position=0 value=11111100:00000000:00000000:00111111", ss.toString());
    assertFalse(ss.isValid());

    assertEquals("11110000:00000000:00000000:00111100", Utils.toBitsSep(ss.next()));
    assertEquals("position=1 value=11110000:00000000:00000000:00111100", ss.toString());
    assertFalse(ss.isValid());

    assertEquals("11000000:00000000:00000000:00110001", Utils.toBitsSep(ss.next()));
    assertEquals("position=2 value=11000000:00000000:00000000:00110001", ss.toString());
    assertFalse(ss.isValid());

    assertEquals("00000000:00000000:00000000:00000110", Utils.toBitsSep(ss.next()));
    assertEquals("position=3 value=00000000:00000000:00000000:00000110", ss.toString());
    assertTrue(ss.isValid());

    assertEquals("11111100:00000000:00000000:00111111", Utils.toBitsSep(ss.next()));
    assertEquals("position=4 value=11111100:00000000:00000000:00111111", ss.toString());
    assertFalse(ss.isValid());

    assertEquals("11110000:00000000:00000000:00111100", Utils.toBitsSep(ss.next()));
    assertEquals("position=5 value=11110000:00000000:00000000:00111100", ss.toString());
    assertFalse(ss.isValid());

    assertEquals("11000000:00000000:00000000:00110001", Utils.toBitsSep(ss.next()));
    assertEquals("position=6 value=11000000:00000000:00000000:00110001", ss.toString());
    assertFalse(ss.isValid());

    assertEquals("00000000:00000000:00000000:00000110", Utils.toBitsSep(ss.next()));
    assertEquals("position=7 value=00000000:00000000:00000000:00000110", ss.toString());
    assertTrue(ss.isValid());

    assertEquals("11111100:00000000:00000000:00111111", Utils.toBitsSep(ss.next()));
    assertEquals("position=8 value=11111100:00000000:00000000:00111111", ss.toString());
    assertFalse(ss.isValid());

    assertEquals("11111100:00000000:00000000:00111111", Utils.toBitsSep(ss.next()));
    assertEquals("position=9 value=11111100:00000000:00000000:00111111", ss.toString());
    assertFalse(ss.isValid());
  }
}
