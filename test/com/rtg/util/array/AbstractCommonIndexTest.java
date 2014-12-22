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

import static com.rtg.util.StringUtils.LS;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractCommonIndexTest extends TestCase {

  protected final void checkSwap(final AbstractIndex index) {
    assertTrue(index.globalIntegrity());
    for (long i = 0; i < index.length(); i++) {
      assertEquals(0, index.get(i));
    }
    // test toString() when all zeroes
    assertEquals("Index [10]" + LS, index.toString());

    index.set(3, 6);
    assertEquals(6, index.get(3));
    index.set(9, 5);
    index.swap(9, 3);
    assertEquals(0, index.get(2));
    assertEquals(5, index.get(3));
    assertEquals(0, index.get(4));
    assertEquals(0, index.get(8));
    assertEquals(6, index.get(9));

    // test toString()
    String str = index.toString().replaceAll("  *", " ");
    assertEquals("Index [10]" + LS + "[0] 0, 0, 0, 5, 0, 0, 0, 0, 0, 6" + LS, str);
    assertTrue(index.globalIntegrity());

    // test toString with a range
    StringBuilder buf = new StringBuilder();
    index.toString(buf, 7, index.length());
    str = buf.toString().replaceAll("  *", " ");
    assertEquals("[7] 0, 0, 6" + LS, str);

    // test dumpString with a range
    buf = new StringBuilder();
    index.dumpString(buf, 7, index.length());
    str = buf.toString().replaceAll("  *", " ");
    assertEquals("[7]"
        + " 00000000:00000000:00000000:00000000:00000000:00000000:00000000:00000000"
        + " 00000000:00000000:00000000:00000000:00000000:00000000:00000000:00000000"
        + " 00000000:00000000:00000000:00000000:00000000:00000000:00000000:00000110" + LS, str);
  }
}
