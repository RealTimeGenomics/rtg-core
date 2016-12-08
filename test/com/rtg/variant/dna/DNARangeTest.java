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
public class DNARangeTest extends IntegerRangeTest {

  @Override
  protected IntegerRange getRange() {
    return new DNARange(0, 4);
  }

  public void testGeneral() {
    final IntegerRange ir = getRange();
    if (ir.hasInvalid()) {
      assertEquals(-1, ir.invalid());
    }
    assertEquals(0, ir.low());
  }

  public void testToChar() {
    final DNARange ir = (DNARange) getRange();
    for (int i = ir.low(); i <= ir.high(); ++i) {
      assertEquals(ir.toString(i), "" + ir.toChar(i));
      assertEquals(i, ir.valueOf(ir.toChar(i)));
    }
  }

  public void testValidExplicit() {
    assertTrue(DNARange.valid(0, DNARange.N, DNARange.T));
    assertTrue(DNARange.valid(1, DNARange.N, DNARange.T));
    assertTrue(DNARange.valid(1, DNARange.N, DNARange.T));

    assertFalse(DNARange.valid(-1, DNARange.N, DNARange.T));
    assertFalse(DNARange.valid(6, DNARange.N, DNARange.T));
  }

  public void testComplement() {
    assertEquals(DNARange.N, DNARange.complement(DNARange.N));
    assertEquals(DNARange.T, DNARange.complement(DNARange.A));
    assertEquals(DNARange.G, DNARange.complement(DNARange.C));
    assertEquals(DNARange.C, DNARange.complement(DNARange.G));
    assertEquals(DNARange.A, DNARange.complement(DNARange.T));
  }
}
