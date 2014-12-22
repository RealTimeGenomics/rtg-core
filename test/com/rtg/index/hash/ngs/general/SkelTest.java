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
