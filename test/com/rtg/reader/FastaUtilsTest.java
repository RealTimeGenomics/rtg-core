/*
 * Copyright (c) 2015. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.reader;

import java.util.Arrays;

import junit.framework.TestCase;

/**
 * JUnit tests for the ReaderUtils class.
 */
public class FastaUtilsTest extends TestCase {

  public void testConversion() {
    assertNull(FastaUtils.asciiToRawQuality((String) null));
    assertNull(FastaUtils.asciiToRawQuality((char[]) null));
    assertNull(FastaUtils.rawToAsciiQuality(null));

    final String qualities = "!\"#$";
    assertTrue(Arrays.equals(new byte[]{0, 1, 2, 3}, FastaUtils.asciiToRawQuality(qualities)));
    assertTrue(Arrays.equals(new byte[]{0, 1, 2, 3}, FastaUtils.asciiToRawQuality(qualities.toCharArray())));

    assertEquals(qualities, new String(FastaUtils.rawToAsciiQuality(FastaUtils.asciiToRawQuality(qualities))));
  }
}
