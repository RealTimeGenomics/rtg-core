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

import junit.framework.TestCase;

/**
 */
public class VariantLocusTest extends TestCase {

  public void testEnough() {
    final String refName = "blkajd";
    final int start = 394;
    final int end = 398;
    final String ref = "AAAT";
    final char prevRefNt = '1';
    final VariantLocus locus = new VariantLocus(refName, start, end, ref, prevRefNt);

    assertEquals(refName, locus.getSequenceName());
    assertEquals(start, locus.getStart());
    assertEquals(end, locus.getEnd());
    assertEquals(ref, locus.getRefNts());
    assertEquals(prevRefNt, locus.getPreviousRefNt());

    assertEquals("start=" + start + " end=" + end, locus.toString());
  }
}
