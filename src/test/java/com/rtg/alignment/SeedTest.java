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
public class SeedTest extends TestCase {

  public void test() {
    final Seed sd = new Seed(4);
    sd.integrity();
    assertEquals("size=4 mask=11111111:00000000:00000000:11111111", sd.toString());
    assertEquals("11111111:00000000:00000000:11111111", Utils.toBitsSep(sd.init()));
    assertEquals(4, sd.size());

    assertEquals(0, sd.next(0, (byte) 1));
    assertEquals(1, sd.next(0, (byte) 2));
    assertEquals(2, sd.next(0, (byte) 3));
    assertEquals(3, sd.next(0, (byte) 4));
    assertEquals("11111111:00000000:00000000:11111111", Utils.toBitsSep(sd.next(0, (byte) 0)));

    int seed = 0;
    for (byte b = 1; b <= 4; ++b) {
      seed = sd.next(seed, b);
      assertTrue(sd.isValid(seed));
    }
    assertEquals("00000000:00000000:00000000:00011011", Utils.toBitsSep(seed));

    seed = sd.next(seed, (byte) 0);
    assertFalse(sd.isValid(seed));
    assertEquals("11111111:00000000:00000000:11111111", Utils.toBitsSep(seed));

    seed = sd.next(seed, (byte) 3);
    assertFalse(sd.isValid(seed));
    assertEquals("11111100:00000000:00000000:11111110", Utils.toBitsSep(seed));

    seed = sd.next(seed, (byte) 3);
    assertFalse(sd.isValid(seed));
    assertEquals("11110000:00000000:00000000:11111010", Utils.toBitsSep(seed));

    seed = sd.next(seed, (byte) 3);
    assertFalse(sd.isValid(seed));
    assertEquals("11000000:00000000:00000000:11101010", Utils.toBitsSep(seed));

    seed = sd.next(seed, (byte) 3);
    assertTrue(sd.isValid(seed));
    assertEquals("00000000:00000000:00000000:10101010", Utils.toBitsSep(seed));
  }

  public void testBad() {
    try {
      new Seed(8);
      fail();
    } catch (final Exception e) {
      //expected
    }
  }
}
