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
package com.rtg.util.array;


import junit.framework.TestCase;

/**
 */
public class ArrayHandleTest extends TestCase {

  /**
   * Test method for {@link com.rtg.util.array.ArrayHandle}.
   */
  public final void test() {
    final ArrayHandle ah = new ArrayHandle(ArrayType.LONG, 100);
    assertEquals(ArrayType.LONG, ah.type());
    assertEquals(100, ah.length());
    final CommonIndex is = ah.createUnsigned();
    assertEquals(100, is.length());
    assertEquals(is.bytes(), ah.bytes());
    final CommonIndex us = ah.createUnsigned();
    assertEquals(100, us.length());
    assertEquals(us.bytes(), ah.bytes());
  }
}

