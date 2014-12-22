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
public class FloatValueTest extends TestCase {

  /**
   */
  public FloatValueTest(final String name) {
    super(name);
  }

  public static Test suite() {
    final TestSuite suite = new TestSuite();
    suite.addTest(new TestSuite(FloatValueTest.class));
    return suite;
  }

  public void test() {
    FloatValue val = new FloatValue(5.23f);

    assertEquals("5.23", val.toString());
    StringBuilder sb = new StringBuilder();
    val.toString(sb);
    assertEquals("5.23", sb.toString());

    assertEquals(Integer.MAX_VALUE, val.maxLength());
    assertEquals(Float.floatToRawIntBits(5.23f), val.hashCode());

    assertTrue(val.equals(val));
    assertTrue(val.equals(new FloatValue(5.23f)));
    assertFalse(val.equals(null));
    assertEquals(0, val.compareTo(val));
    assertEquals(0, val.compareTo(new FloatValue(5.23f)));
    assertEquals(-1, val.compareTo(new FloatValue(5.24f)));
    assertEquals(1, val.compareTo(new FloatValue(5.22f)));
    assertEquals(1, val.compareTo(NullValue.NULL));
  }
}
