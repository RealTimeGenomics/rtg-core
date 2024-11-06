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
package com.rtg.index.hash.ngs.general;

import com.rtg.util.integrity.Exam;

import junit.framework.TestCase;

/**
 */
public class SkelTest extends TestCase {

  public void test1() {
    final Skel sk = new Skel(7, 4, 6);
    sk.integrity();
    assertEquals("[7...4]<<3", sk.toString());
    assertEquals(4, sk.length());
    assertEquals(7, sk.position());
    assertEquals(6, sk.finalPosition());
    assertEquals(0x78, sk.mask(-1L));
    assertEquals(0x78, sk.mask(-1L, 0));
    assertEquals(0x78, sk.mask(-1L, 1));
    assertEquals(0x78, sk.mask(-1L, -1));
  }

  public void test2() {
    final Skel sk = new Skel(63, 64, 63);
    assertEquals(63, sk.position());
    sk.integrity();
    assertEquals("[63...0]<<0", sk.toString());
    assertEquals(63, sk.finalPosition());
    assertEquals(-1, sk.mask(-1L));
    assertEquals(-1, sk.mask(-1L, 0));
  }

  public void test3() {
    final Skel sk = new Skel(0, 1, 0);
    sk.integrity();
    assertEquals("[0...0]<<0", sk.toString());
    assertEquals(0, sk.finalPosition());
    assertEquals(1, sk.length());
    assertEquals(1, sk.mask(-1L));
  }

  public void test4() {
    final Skel sk = new Skel(63, 1, 63);
    sk.integrity();
    assertEquals("[63...63]<<63", sk.toString());
    assertEquals(1, sk.length());
    assertEquals(Long.MIN_VALUE, sk.mask(-1L));
    assertEquals(Long.MIN_VALUE, sk.mask(-1L, 0));
  }

  //triggered bug in initial testing
  public void test5() {
    final Skel sk = new Skel(7, 4, 6);
    sk.integrity();
    assertEquals("[7...4]<<3", sk.toString());
    assertEquals(4, sk.length());
    assertEquals(0x50, sk.mask(0xAF));
    assertEquals(0x50, sk.mask(0xAF, 0));
    assertEquals(0x28, sk.mask(0xAF, 1));
    assertEquals(0x28, sk.mask(0xAF, -1));
  }

  public void test6() {
    final Skel sk = new Skel(7, 4, 6);
    sk.integrity();
    assertEquals("[7...4]<<3", sk.toString());
    assertEquals(4, sk.length());
    assertEquals(0x78, sk.mask(0xF0));
    assertEquals(0x78, sk.mask(0xF0, 0));
    assertEquals(0x38, sk.mask(0xF0, 1));
    assertEquals(0x70, sk.mask(0xF0, -1));
  }

  public void test7() {
    final Skel sk = new Skel(7, 4, 6);
    sk.integrity();
    assertEquals("[7...4]<<3", sk.toString());
    assertEquals(4, sk.length());
    assertEquals(0x58, sk.mask(0xB0));
    assertEquals(0x58, sk.mask(0xB0, 0));
    assertEquals(0x28, sk.mask(0xB0, 1));
    assertEquals(0x30, sk.mask(0xB0, -1));
  }

  public void testInvalid() {
    checkInvalid(-1, 1, 0);
    checkInvalid(65, 1, 0);
    checkInvalid(2, 4, 0);
    checkInvalid(2, 0, 0);
    checkInvalid(2, 1, -1);
    checkInvalid(2, 1, 3);
    checkInvalid(64, 64, 1);
  }

  private void checkInvalid(final int p, final int ln, final int rt) {
    try {
      new Skel(p, ln, rt).integrity();
      fail();
    } catch (final Exam.ExamException e) {
      //expected
    }
  }
}
