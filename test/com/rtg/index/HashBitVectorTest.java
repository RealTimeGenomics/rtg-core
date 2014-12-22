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
public class HashBitVectorTest extends TestCase {

  public HashBitVectorTest(final String name) {
    super(name);
  }

  public void test0() {
    final HashBitVector bitv = new HashBitVector(1, 0);
    bitv.integrity();
    assertEquals(1, bitv.length());
    assertEquals(4, bitv.bytes());

    assertTrue(!bitv.get(0));
    assertTrue(!bitv.get(1));
    bitv.set(0);
    assertTrue(bitv.get(0));
    assertTrue(bitv.get(1));
    bitv.set(1);
    assertTrue(bitv.get(0));
    assertTrue(bitv.get(1));
    final String expected = ""
      + "HashBitVector bits=1 vectorBits=0 shift=1" + StringUtils.LS
      + "BitVector[1]" + StringUtils.LS
      + "[0]\t 1" + StringUtils.LS
      + StringUtils.LS;
    assertEquals(expected, bitv.toString());
  }


  public void test42() {
    final HashBitVector bitv = new HashBitVector(4, 2);
    bitv.integrity();
    assertTrue(!bitv.get(1));
    assertTrue(!bitv.get(2));
    assertTrue(!bitv.get(15));
    assertTrue(!bitv.get(13));

    bitv.set(0);
    assertTrue(bitv.get(1));
    assertTrue(bitv.get(2));
    assertTrue(!bitv.get(15));
    assertTrue(!bitv.get(13));

    bitv.set(12);
    assertTrue(bitv.get(1));
    assertTrue(bitv.get(2));
    assertTrue(bitv.get(15));
    assertTrue(bitv.get(13));


    final String expected = ""
        + "HashBitVector bits=4 vectorBits=2 shift=2" + StringUtils.LS
        + "BitVector[4]" + StringUtils.LS
        + "[0]\t 1001" + StringUtils.LS
        + StringUtils.LS;
    assertEquals(expected, bitv.toString());
  }

  public void testDirect() {
    final HashBitVector bitv = new HashBitVector(4, 2);
    bitv.integrity();
    assertTrue(!bitv.getDirect(0));
    assertTrue(!bitv.getDirect(1));
    assertTrue(!bitv.getDirect(2));
    assertTrue(!bitv.getDirect(3));

    bitv.setDirect(0);
    assertTrue(bitv.getDirect(0));
    assertTrue(!bitv.getDirect(1));
    assertTrue(!bitv.getDirect(2));
    assertTrue(!bitv.getDirect(3));

    final String expected = ""
        + "HashBitVector bits=4 vectorBits=2 shift=2" + StringUtils.LS
        + "BitVector[4]" + StringUtils.LS
        + "[0]\t 1000" + StringUtils.LS
        + StringUtils.LS;
    assertEquals(expected, bitv.toString());

    bitv.resetDirect(0);
    assertTrue(!bitv.getDirect(0));
    assertTrue(!bitv.getDirect(1));
    assertTrue(!bitv.getDirect(2));
    assertTrue(!bitv.getDirect(3));

  }

  public void test22() {
    final HashBitVector bitv = new HashBitVector(2, 2);
    bitv.integrity();
    assertTrue(!bitv.get(1));
    assertTrue(!bitv.get(2));

    bitv.set(0);
    assertTrue(bitv.get(0));
    assertTrue(!bitv.get(2));

    bitv.set(2);
    assertTrue(bitv.get(0));
    assertTrue(bitv.get(2));


    final String expected = ""
      + "HashBitVector bits=2 vectorBits=2 shift=0" + StringUtils.LS
      + "BitVector[4]" + StringUtils.LS
      + "[0]\t 1010" + StringUtils.LS
      + StringUtils.LS;
    assertEquals(expected, bitv.toString());
  }

  public void test64x2() {
    final HashBitVector bitv = new HashBitVector(64, 2);
    bitv.integrity();
    assertTrue(!bitv.get(1));
    assertTrue(!bitv.get(-1L));

    bitv.set(0);
    assertTrue(bitv.get(0));
    assertTrue(!bitv.get(-1L));

    bitv.set(-1L);
    assertTrue(bitv.get(0));
    assertTrue(bitv.get(-1L));

    final String expected = ""
      + "HashBitVector bits=64 vectorBits=2 shift=62" + StringUtils.LS
      + "BitVector[4]" + StringUtils.LS
      + "[0]\t 1001" + StringUtils.LS
      + StringUtils.LS;
    assertEquals(expected, bitv.toString());
  }

  public void test1x4() {
    //bits > vectorbits
    final HashBitVector bitv = new HashBitVector(1, 4);
    bitv.integrity();
    assertTrue(!bitv.get(1));
    assertTrue(!bitv.get(0));

    bitv.set(0);
    assertTrue(bitv.get(0));
    assertTrue(!bitv.get(1));

    bitv.set(1);
    assertTrue(bitv.get(0));
    assertTrue(bitv.get(1));

    final String expected = ""
      + "HashBitVector bits=1 vectorBits=4 shift=0" + StringUtils.LS
      + "BitVector[16]" + StringUtils.LS
      + "[0]\t 1100000000 000000" + StringUtils.LS
      + StringUtils.LS;
    assertEquals(expected, bitv.toString());
  }

  public void testExceptions() {
    try {
      new HashBitVector(1, -1);
      fail("Exception expected");
    } catch (final Exception e) {
      assertEquals("Invalid bit parameters bits=1 vectorBits=-1", e.getMessage());
    }
    try {
      new HashBitVector(-1, 1);
      fail("Exception expected");
    } catch (final Exception e) {
      assertEquals("Invalid bit parameters bits=-1 vectorBits=1", e.getMessage());
    }
    try {
      new HashBitVector(65, 3);
      fail("Exception expected");
    } catch (final Exception e) {
      assertEquals("Invalid bit parameters bits=65 vectorBits=3", e.getMessage());
    }

  }
}

