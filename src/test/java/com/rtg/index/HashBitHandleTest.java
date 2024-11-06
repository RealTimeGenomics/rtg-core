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
package com.rtg.index;

import com.rtg.util.integrity.Exam;

import junit.framework.TestCase;

/**
 */
public class HashBitHandleTest extends TestCase {

  public void test1() {
    final HashBitHandle hbh = new HashBitHandle(5, 7);
    assertEquals(128, hbh.length());
    assertEquals(7, hbh.bits());
    assertEquals(16, hbh.bytes());
    assertEquals("HashBit Handle hash bits=5 length bits=7", hbh.toString());
    final HashBitVector hbv = hbh.create();
    assertTrue(hbv.toString().startsWith("HashBitVector bits=5 vectorBits=7 shift=0"));
    assertEquals(128, hbv.length());
    assertEquals(16, hbv.bytes());
  }

  public void test2() {
    final HashBitHandle hbh = new HashBitHandle(8, 7);
    assertEquals(128, hbh.length());
    assertEquals(7, hbh.bits());
    assertEquals(16, hbh.bytes());
    assertEquals("HashBit Handle hash bits=8 length bits=7", hbh.toString());
    final HashBitVector hbv = hbh.create();
    assertTrue(hbv.toString().startsWith("HashBitVector bits=8 vectorBits=7 shift=1"));
    assertEquals(128, hbv.length());
    assertEquals(16, hbv.bytes());
  }

  public void test3() {
    final HashBitHandle hbh = new HashBitHandle(64, 62);
    assertEquals(1L << 62, hbh.length());
    assertEquals(62, hbh.bits());
    assertEquals(1L << 59, hbh.bytes());
    assertEquals("HashBit Handle hash bits=64 length bits=62", hbh.toString());
  }

  public void test4() {
    final HashBitHandle hbh = new HashBitHandle(62, 62);
    assertEquals(1L << 62, hbh.length());
    assertEquals(62, hbh.bits());
    assertEquals(1L << 59, hbh.bytes());
    assertEquals("HashBit Handle hash bits=62 length bits=62", hbh.toString());
  }

  public void test5() {
    final HashBitHandle hbh = new HashBitHandle(61, 62);
    assertEquals(1L << 62, hbh.length());
    assertEquals(62, hbh.bits());
    assertEquals(1L << 59, hbh.bytes());
    assertEquals("HashBit Handle hash bits=61 length bits=62", hbh.toString());
  }

  public void testBad() {
    try {
      new HashBitHandle(64, 63);
      fail();
    } catch (final Exam.ExamException e) {
      // expected
    }
    try {
      new HashBitHandle(63, 63);
      fail();
    } catch (final Exam.ExamException e) {
      // expected
    }

    try {
      new HashBitHandle(-1, 1);
      fail();
    } catch (final Exam.ExamException e) {
      // expected
    }
    new HashBitHandle(1, 0);

  }
}

