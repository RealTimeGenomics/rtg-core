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
package com.rtg.variant.match;


import junit.framework.TestCase;

/**
 * Test class
 */
public class AlleleAsReadMatchTest extends TestCase {

  public AlleleAsReadMatchTest(final String testName) {
    super(testName);
  }

  private AlleleAsReadMatch mMatch = null;
  @Override
  public void setUp() {
    mMatch = new AlleleAsReadMatch(new byte[] {4, 3, 2, 1, 0}, 20.0);
  }

  @Override
  public void tearDown() {
    mMatch = null;
  }

  public void testIsFixedLeft() {
    assertTrue(mMatch.isFixedLeft());
  }

  public void testIsFixedRight() {
    assertTrue(mMatch.isFixedRight());
  }

  public void testLength() {
    assertEquals(5, mMatch.length());
  }

  public void testMapQuality() {
    assertEquals(0.0, mMatch.mapError());
  }

  public void testQuality() {
    assertEquals(20.0, mMatch.baseError(3));
  }

  public void testRead() {
    assertEquals(3, mMatch.read(0));
    assertEquals(2, mMatch.read(1));
    assertEquals(1, mMatch.read(2));
    assertEquals(0, mMatch.read(3));
  }

  public void testReadString() {
    assertEquals("TGCAN", mMatch.readString());
  }

  public void testEmpty() {
    final Match ref = new AlleleAsReadMatch(new byte[] {}, 20.0);
    assertEquals("", ref.readString());
  }
}
