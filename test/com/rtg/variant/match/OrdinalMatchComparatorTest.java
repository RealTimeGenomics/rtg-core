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
package com.rtg.variant.match;


import junit.framework.TestCase;

/**
 */
public class OrdinalMatchComparatorTest extends TestCase {

  public void test() {
    final OrdinalMatchComparator comparator = new OrdinalMatchComparator();
    final Match snp1 = new AlignmentMatch(null, "AC", "!!", 0, 0, 2, 0);
    final Match snp2 = new AlignmentMatch(null, "CC", "!!", 0, 0, 2, 0);
    assertEquals(0, comparator.compare(snp1, snp1));
    assertEquals(-2, comparator.compare(snp1, snp2));
    assertEquals(2, comparator.compare(snp2, snp1));
  }

}
