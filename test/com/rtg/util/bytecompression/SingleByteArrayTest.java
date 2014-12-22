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
package com.rtg.util.bytecompression;

import java.io.ByteArrayInputStream;

import junit.framework.TestCase;

/**
 * Tests for the corresponding class
 *
 */

public class SingleByteArrayTest extends TestCase {

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(SingleByteArrayTest.class);
  }

  public void test() throws Exception {
    final SingleByteArray sba = new SingleByteArray(10);
    assertEquals(10, sba.length());
    sba.set(2, new byte[]{3}, 1);
    assertEquals(3, MultiByteArrayTest.get1(sba, 2));

    final byte[] b = new byte[3];
    sba.get(b, 1, 2);
    assertEquals(0, b[0]);
    assertEquals(3, b[1]);
    assertEquals(0, b[2]);
    assertEquals(0, sba.get(0));
    assertEquals(0, sba.get(1));
    assertEquals(3, sba.get(2));
    try {
      sba.load(new ByteArrayInputStream(b), 0, 55);
      fail("expected index out of bounds exception");
    } catch (final IndexOutOfBoundsException ioe) {
      //expected
    }
    sba.load(new ByteArrayInputStream(b), 0, 2);
  }

}
