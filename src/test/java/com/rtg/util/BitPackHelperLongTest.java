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
package com.rtg.util;

import junit.framework.TestCase;

/**
 */
public class BitPackHelperLongTest extends TestCase {

  public void test() {
    try {
      new BitPackHelperLong(32, 16, 16, 1);
      fail();
    } catch (final IllegalArgumentException e) {

    }
    final BitPackHelperLong bit = new BitPackHelperLong(14, 15, 16, 19);
    long value = bit.packValues(new long[] {34, 53, 65535, 524287});
    assertEquals(34, bit.getField(0, value));
    assertEquals(53, bit.getField(1, value));
    assertEquals(65535, bit.getField(2, value));
    assertEquals(524287, bit.getField(3, value));
    value = bit.packValues(new long[] {12, 12, 12, 12});
    assertEquals(12, bit.getField(0, value));
    assertEquals(12, bit.getField(1, value));
    assertEquals(12, bit.getField(2, value));
    assertEquals(12, bit.getField(3, value));
    value = bit.packValues(new long[] {205, 1105, 405, 905});
    assertEquals(205, bit.getField(0, value));
    assertEquals(1105, bit.getField(1, value));
    assertEquals(405, bit.getField(2, value));
    assertEquals(905, bit.getField(3, value));
  }
}
