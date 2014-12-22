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
package com.rtg.mode;

import static com.rtg.mode.SequenceType.DNA;
import static com.rtg.mode.SequenceType.PROTEIN;

import java.io.ObjectStreamException;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class SequenceTypeTest extends TestCase {

  /**
   * Test method for {@link com.rtg.mode.SequenceType}.
   * @throws ObjectStreamException
   */
  public final void test() throws ObjectStreamException {
    TestUtils.testPseudoEnum(SequenceType.class, "[DNA, PROTEIN]");
    assertEquals(SequenceType.DNA, SequenceType.DNA.readResolve());
    assertEquals(SequenceType.PROTEIN, SequenceType.PROTEIN.readResolve());
  }

  /**
   * Test method for {@link com.rtg.mode.SequenceType#numberKnownCodes()}.
   */
  public final void testNumberKnownCodes() {
    assertEquals(4, DNA.numberKnownCodes());
    assertEquals(5, DNA.numberCodes());
    assertEquals(com.rtg.mode.DNA.values().length, DNA.firstValid() + DNA.numberKnownCodes());
    assertEquals(20, PROTEIN.numberKnownCodes());
    assertEquals(22, PROTEIN.numberCodes());
    assertEquals(com.rtg.mode.Protein.values().length, PROTEIN.firstValid() + PROTEIN.numberKnownCodes());
    assertEquals("DNA", SequenceType.DNA.name());
  }

  /**
   * Test method for {@link com.rtg.mode.SequenceType#bits()}.
   */
  public final void testBits() {
    check(DNA);
    check(PROTEIN);
  }

  private void check(final SequenceType st) {
    final int b = st.bits();
    assertTrue(b < 64);
    final long cb = 1 << b;
    final int n = st.numberKnownCodes();
    assertTrue(cb >= n);
    final int bm = b - 1;
    final long cbm = 1 << bm;
    assertTrue(cbm < n);
  }
}

