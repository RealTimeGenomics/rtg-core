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

package com.rtg.variant.avr;

import junit.framework.TestCase;

/**
 */
public class IntegerConverterTest extends TestCase {
  public void testConvert() {
    final IntegerConverter converter = new IntegerConverter();
    assertNull(converter.stringToObjectOfType(null));
     try {
      converter.stringToObjectOfType("apple");
      fail("converted apple");
    } catch (NumberFormatException nfe) {
       // expected
    }
    assertEquals(1234, converter.stringToObjectOfType("1234"));
    assertTrue(converter.stringToObjectOfType("1234") instanceof Integer);
    try {
      converter.stringToObjectOfType("12.34");
      fail("converted 12.34");
    } catch (NumberFormatException nfe) {
      // expected
    }
  }

}
