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
package com.rtg.util.array.bitindex;

import static com.rtg.util.StringUtils.LS;

import junit.framework.TestCase;


/**
 */
public class BitCreateTest extends TestCase {

  public void testTooSmall() {
    try {
      BitCreate.createIndex(10L, 0);
      fail("exception expected");
    } catch (final Exception e) {
      assertEquals("Illegal bits value=" + 0, e.getMessage());
    }
  }

  public void testTooBig() {
    try {
      BitCreate.createIndex(10L, 65);
      fail("exception expected");
    } catch (final Exception e) {
      assertEquals("Illegal bits value=" + 65, e.getMessage());
    }
  }

  public void testOK() {
    final BitIndex index = BitCreate.createIndex(10L, 64);
    assertEquals("Index [10]" + LS, index.toString());
  }

}
