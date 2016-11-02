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

package com.rtg.variant;

import java.util.Arrays;

import junit.framework.TestCase;

/**
 */
public class InterpolateArrayTest extends TestCase {

  public void testSimple() {
    final int[] curve = {-1, 4, -1, 6, -1, 10, -1};
    final int[] expected = {-1, 4, 5, 6, 8, 10, -1};


    new InterpolateArray(curve).interpolate();
    final String exp = Arrays.toString(expected);
    final String actual = Arrays.toString(curve);
    assertTrue(String.format("Expected <%s> but was <%s>", exp, actual), Arrays.equals(expected, curve));
  }

  public void testFill() {
    final int[] curve = {-1, 4, -1, 6, -1, 10, -1};
    final int[] defaults = {99, 98, 97, 96, 95, 94, 93};

    final int[] expected = {99, 4, 97, 6, 95, 10, 93};

    new InterpolateArray(curve).fill(defaults);

    final String exp = Arrays.toString(expected);
    final String actual = Arrays.toString(curve);
    assertTrue(String.format("Expected <%s> but was <%s>", exp, actual), Arrays.equals(expected, curve));

  }

}
