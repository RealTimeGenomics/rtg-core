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
public class BitPack2IntoLongTest extends TestCase {

  public void test() {
    try {
      new BitPack2IntoLong(64, 1);
      fail();
    } catch (final IllegalArgumentException e) {
    }
    try {
      new BitPack2IntoLong(60, 0);
      fail();
    } catch (final IllegalArgumentException e) {
    }
    try {
      new BitPack2IntoLong(0, 60);
      fail();
    } catch (final IllegalArgumentException e) {
    }
    final BitPack2IntoLong bit = new BitPack2IntoLong(13, 8);
    long value = bit.packValues(34, 90);
    assertEquals(34, bit.getField(0, value));
    assertEquals(90, bit.getField(1, value));
    value = bit.packValues(12, 12);
    assertEquals(12, bit.getField(0, value));
    assertEquals(12, bit.getField(1, value));
    value = bit.packValues(25, 11);
    assertEquals(25, bit.getField(0, value));
    assertEquals(11, bit.getField(1, value));
  }
}
