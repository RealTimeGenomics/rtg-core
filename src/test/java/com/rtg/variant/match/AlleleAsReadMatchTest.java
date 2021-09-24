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
