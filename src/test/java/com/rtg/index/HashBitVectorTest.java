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

