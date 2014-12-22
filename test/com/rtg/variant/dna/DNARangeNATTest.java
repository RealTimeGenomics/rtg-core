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

package com.rtg.variant.dna;

import com.rtg.util.integrity.IntegerRange;
import com.rtg.util.integrity.IntegerRangeTest;

/**
 */
public class DNARangeNATTest extends IntegerRangeTest {

  @Override
  protected IntegerRange getRange() {
    return DNARangeNAT.DNA;
  }

  public void testLocal() {
    DNARangeNAT.DNA.valid(DNARangeNAT.A);
    DNARangeNAT.DNA.valid(DNARangeNAT.C);
    DNARangeNAT.DNA.valid(DNARangeNAT.G);
    DNARangeNAT.DNA.valid(DNARangeNAT.T);
    assertTrue(DNARangeNAT.DNA.hasInvalid());
    assertEquals(0, DNARangeNAT.DNA.low());
    assertEquals(3, DNARangeNAT.DNA.high());
  }

}
