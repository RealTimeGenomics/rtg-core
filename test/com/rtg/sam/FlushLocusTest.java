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

package com.rtg.sam;

import junit.framework.TestCase;

/**
 */
public class FlushLocusTest extends TestCase {

  public void test() {
    final FlushLocus fla = new FlushLocus(0, 10);
    assertEquals("[0,10)", fla.toString());
  }

  public void testEqHash() {
    final FlushLocus fla = new FlushLocus(0, 10);
    final FlushLocus flb = new FlushLocus(5, 15);

    assertEquals(202, fla.hashCode());
    assertEquals(fla, fla);
    assertFalse(fla.equals(flb));
    assertFalse(flb.equals(fla));

    final FlushLocus flx = new FlushLocus(0, 11);
    final FlushLocus fly = new FlushLocus(1, 10);
    assertFalse(fla.equals(flx));
    assertFalse(flx.equals(fla));
    assertFalse(fla.equals(fly));
    assertFalse(fly.equals(fla));

  }

  public void testJoinable() {
    checkJoinable(new FlushLocus(0, 10), new FlushLocus(0, 10), new FlushLocus(0, 10));
    checkJoinable(new FlushLocus(0, 10), new FlushLocus(5, 15), new FlushLocus(0, 15));
    checkJoinable(new FlushLocus(0, 10), new FlushLocus(10, 20), new FlushLocus(0, 20));
    checkJoinable(new FlushLocus(0, 10), new FlushLocus(15, 25), null);

    assertEquals(202, new FlushLocus(0, 10).hashCode());
    assertEquals(new FlushLocus(0, 10), new FlushLocus(0, 10));
    assertFalse(new FlushLocus(0, 10).equals(new FlushLocus(5, 15)));
    assertFalse(new FlushLocus(5, 15).equals(new FlushLocus(0, 10)));

    final FlushLocus flx = new FlushLocus(0, 11);
    final FlushLocus fly = new FlushLocus(1, 10);
    assertFalse(new FlushLocus(0, 10).equals(flx));
    assertFalse(flx.equals(new FlushLocus(0, 10)));
    assertFalse(new FlushLocus(0, 10).equals(fly));
    assertFalse(fly.equals(new FlushLocus(0, 10)));

  }

  private void checkJoinable(final FlushLocus fla, final FlushLocus flb, final FlushLocus flc) {
    final boolean exp = flc != null;
    assertEquals(exp, fla.isJoinable(flb));
    assertEquals(exp, flb.isJoinable(fla));
    if (exp) {
      fla.join(flb);
      assertEquals(fla.mStart, flc.mStart);
      assertEquals(fla.mEnd, flc.mEnd);
    }
  }
}
