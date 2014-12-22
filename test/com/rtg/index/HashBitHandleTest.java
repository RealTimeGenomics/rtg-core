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

