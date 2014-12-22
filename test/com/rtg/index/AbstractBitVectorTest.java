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

import com.rtg.util.StringUtils;

import junit.framework.TestCase;


/**
 * JUnit tests for BitVector.
 *
 */
public abstract class AbstractBitVectorTest extends TestCase {

  protected abstract AbstractBitVector getBitVector(long length);

  public void testCreate() {
      final AbstractBitVector bv = AbstractBitVector.createBitVector(10);
      assertTrue(bv instanceof BitVector);
  }

  public void test0() {
    final AbstractBitVector bitv = getBitVector(0);
    bitv.integrity();
    assertEquals(0, bitv.length());
    assertEquals(0, bitv.bytes());
    try {
      bitv.get(0);
      fail("Exception expected");
    } catch (final ArrayIndexOutOfBoundsException e) {
      assertEquals("0:0", e.getMessage());
    }
    try {
      bitv.set(0);
      fail("Exception expected");
    } catch (final ArrayIndexOutOfBoundsException e) {
      assertEquals("0:0", e.getMessage());
    }
    try {
      bitv.reset(0);
      fail("Exception expected");
    } catch (final ArrayIndexOutOfBoundsException e) {
      assertEquals("0:0", e.getMessage());
    }
    final String expected = ""
      + "BitVector[0]" + StringUtils.LS
      + StringUtils.LS;
    assertEquals(expected, bitv.toString());
  }


  public void test1() {
    final AbstractBitVector bitv = getBitVector(1);
    bitv.integrity();
    assertEquals(1, bitv.length());
    assertEquals(4, bitv.bytes());
    assertTrue(!bitv.get(0));
    bitv.set(0);
    assertTrue(bitv.get(0));
    final String expected = ""
      + "BitVector[1]" + StringUtils.LS
      + "[0]\t 1" + StringUtils.LS
      + StringUtils.LS;
    assertEquals(expected, bitv.toString());

    bitv.reset(0);
    assertTrue(!bitv.get(0));
  }

  public void testSetAndGet() {
    final AbstractBitVector bitv = getBitVector(100);
    bitv.integrity();
    assertEquals(100, bitv.length());
    assertEquals(16, bitv.bytes());
    for (int i = 0; i < 100; i++) {
      assertTrue(!bitv.get(i));
      if (i % 2 == 0) {
        bitv.set(i);
        assertTrue(bitv.get(i));
      }
    }
    for (int i = 0; i < 100; i++) {
      if (i % 2 == 0) {
        assertTrue(bitv.get(i));
      } else {
        assertTrue(!bitv.get(i));
      }
    }
    final String expected = ""
      + "BitVector[100]" + StringUtils.LS
      + "[0]\t 1010101010 1010101010 1010101010 1010101010 1010101010 1010101010 1010101010 1010101010 1010101010 1010101010" + StringUtils.LS
      + StringUtils.LS;
    assertEquals(expected, bitv.toString());
  }

  public void testResetAndGet() {
    final AbstractBitVector bitv = getBitVector(100);
    bitv.integrity();
    assertEquals(100, bitv.length());
    assertEquals(16, bitv.bytes());
    for (int i = 0; i < 100; i++) {
      assertTrue(!bitv.get(i));
      if (i % 2 == 0) {
        bitv.set(i);
        assertTrue(bitv.get(i));
      }
    }
    for (int i = 0; i < 100; i++) {
      if (i % 4 == 0) {
        bitv.reset(i);
      }
    }
    for (int i = 0; i < 100; i++) {
      if (i % 4 == 2) {
        assertTrue(bitv.get(i));
      } else {
        assertTrue(!bitv.get(i));
      }
    }
    final String expected = ""
      + "BitVector[100]" + StringUtils.LS
      + "[0]\t 0010001000 1000100010 0010001000 1000100010 0010001000 1000100010 0010001000 1000100010 0010001000 1000100010" + StringUtils.LS
      + StringUtils.LS;
    assertEquals(expected, bitv.toString());
  }

  public void testToString123() {
    final AbstractBitVector bitv = getBitVector(123);
    bitv.integrity();
    final String expected = ""
      + "BitVector[123]" + StringUtils.LS
      + "[0]\t 0000000000 0000000000 0000000000 0000000000 0000000000 0000000000 0000000000 0000000000 0000000000 0000000000" + StringUtils.LS
      + "[100]\t 0000000000 0000000000 000" + StringUtils.LS
      + StringUtils.LS;
    assertEquals(expected, bitv.toString());
  }

  public void testToString20() {
    final AbstractBitVector bitv = getBitVector(20);
    bitv.integrity();
    final String expected = ""
      + "BitVector[20]" + StringUtils.LS
      + "[0]\t 0000000000 0000000000" + StringUtils.LS
      + StringUtils.LS;
    assertEquals(expected, bitv.toString());
  }

  public void testOutOfBounds() {
    final AbstractBitVector bitv = getBitVector(100);
    try {
      bitv.get(100);
      fail("Exception expected");
    } catch (final ArrayIndexOutOfBoundsException e) {
      assertEquals("100:100", e.getMessage());
    }
    try {
      bitv.set(100);
      fail("Exception expected");
    } catch (final ArrayIndexOutOfBoundsException e) {
      assertEquals("100:100", e.getMessage());
    }
    try {
      bitv.get(-1);
      fail("Exception expected");
    } catch (final ArrayIndexOutOfBoundsException e) {
      assertEquals("-1:100", e.getMessage());
    }
    try {
      bitv.set(-1);
      fail("Exception expected");
    } catch (final ArrayIndexOutOfBoundsException e) {
      assertEquals("-1:100", e.getMessage());
    }
  }

  public void testExceptions() {
    try {
      getBitVector(-1);
      fail("Exception expected");
    } catch (final Exception e) {
      assertEquals("Negative length:-1", e.getMessage());
    }
    getBitVector(0);
  }
}

