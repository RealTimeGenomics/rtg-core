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

package com.rtg.variant.sv;


/**
 */
public class SamArrayTest extends AbstractSamCountsTest {

  @Override
  protected SamCounts getCounts(int length) {
    return new SamArray(length);
  }

  @Override
  protected void check10Inc(final SamCounts sa) {
    sa.increment(0);
    ((SamArray) sa).increment(4, 2.0);
  }

  public void testReverse() {
    final SamArray sa = new SamArray(5);
    for (int i = 0; i < 5; i++) {
      sa.increment(i, i + 1);
    }
    assertEquals("[" + 5.0 + ", " + 5.0 + ", " + 4.0 + ", " + 3.0 + ", " + 2.0 + "]", sa.reverse(3, 0).toString());
    assertEquals("[" + 5.0 + ", " + 4.0 + ", " + 3.0 + ", " + 2.0 + ", " + 1.0 + "]", sa.reverse(3, 1).toString());
    assertEquals("[" + 4.0 + ", " + 3.0 + ", " + 2.0 + ", " + 1.0 + ", " + 1.0 + "]", sa.reverse(3, 2).toString());
  }
}
