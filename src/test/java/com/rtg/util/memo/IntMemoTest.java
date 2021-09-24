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
package com.rtg.util.memo;

import static com.rtg.util.StringUtils.LS;

import junit.framework.TestCase;

/**
 */
public class IntMemoTest extends TestCase {

  public void test() {
    final IntFunction function = new IntFunction() {
      @Override
      public int fn(int i) {
        return i + 42;
      }
    };
    final IntMemo im = new IntMemo(function);
    im.globalIntegrity();
    assertEquals("IntMemo empty" + LS, im.toString());
    assertEquals(42, im.fn(0));
    im.globalIntegrity();
    assertEquals(44, im.fn(2));
    im.globalIntegrity();
    assertEquals(37, im.fn(-5));
    im.globalIntegrity();
    assertEquals("IntMemo [-5..2]" + LS + "37, 38, 39, 40, 41, 42, 43, 44" + LS, im.toString());

    assertEquals(35, im.fn(-7));
    im.globalIntegrity();
    assertEquals(46, im.fn(4));
    im.globalIntegrity();
    assertEquals("IntMemo [-7..4]" + LS + "35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46" + LS, im.toString());
    assertEquals(42, im.fn(0));
    assertEquals(44, im.fn(2));
    assertEquals(37, im.fn(-5));
  }
}
