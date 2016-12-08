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
public class GeneralDNARangeTest extends IntegerRangeTest {

  @Override
  protected IntegerRange getRange() {
    return new GeneralDNARange("XYZ", false);
  }

  public void testGeneral() {
    final IntegerRange ir = getRange();
    if (ir.hasInvalid()) {
      assertEquals(-1, ir.invalid());
    }
    assertEquals(0, ir.low());
  }

  public void testToChar() {
    final GeneralDNARange ir = (GeneralDNARange) getRange();
    for (int i = ir.low(); i <= ir.high(); ++i) {
      assertEquals(ir.toString(i), "" + ir.toChar(i));
      assertEquals(i, ir.valueOf(ir.toChar(i)));
    }
  }

}
