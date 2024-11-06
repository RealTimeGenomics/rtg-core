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
