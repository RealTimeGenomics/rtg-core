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
package com.rtg.util.format;


import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 */
public class NullValueTest extends TestCase {

  /**
   */
  public NullValueTest(final String name) {
    super(name);
  }

  public static Test suite() {
    final TestSuite suite = new TestSuite();
    suite.addTest(new TestSuite(NullValueTest.class));
    return suite;
  }

  public void test() {
    assertEquals("--", NullValue.NULL.toString());

    StringBuilder sb = new StringBuilder();
    NullValue.NULL.toString(sb);
    assertEquals("--", sb.toString());
    assertEquals(2, NullValue.NULL.maxLength());

    assertEquals(-1, NullValue.NULL.hashCode());

    assertFalse(NullValue.NULL.equals(true));
    assertTrue(NullValue.NULL.equals(NullValue.NULL));

    assertEquals(0, NullValue.NULL.compareTo(NullValue.NULL));
  }
}
