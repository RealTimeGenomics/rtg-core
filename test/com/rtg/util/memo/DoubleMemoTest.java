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
 *
 *
 */
public class DoubleMemoTest extends TestCase {

  public void test() {
    final DoubleFunction function = new DoubleFunction() {
      @Override
      public double fn(int i) {
        return i + 42.3;
      }
    };
    final DoubleMemo im = new DoubleMemo(function);
    im.globalIntegrity();
    assertEquals("IntMemo empty" + LS, im.toString());
    assertEquals(42.3, im.fn(0));
    im.globalIntegrity();
    assertEquals(44.3, im.fn(2));
    im.globalIntegrity();
    assertEquals(37.3, im.fn(-5));
    im.globalIntegrity();
    assertEquals("IntMemo [-5..2]" + LS + "37.3, 38.3, 39.3, 40.3, 41.3, 42.3, 43.3, 44.3" + LS, im.toString());

    assertEquals(35.3, im.fn(-7));
    im.globalIntegrity();
    assertEquals(46.3, im.fn(4));
    im.globalIntegrity();
    assertEquals("IntMemo [-7..4]" + LS + "35.3, 36.3, 37.3, 38.3, 39.3, 40.3, 41.3, 42.3, 43.3, 44.3, 45.3, 46.3" + LS, im.toString());
    assertEquals(42.3, im.fn(0));
    assertEquals(44.3, im.fn(2));
    assertEquals(37.3, im.fn(-5));
  }
}
