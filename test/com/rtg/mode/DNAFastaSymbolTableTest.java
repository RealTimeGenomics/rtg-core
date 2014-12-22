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


import junit.framework.TestCase;

/**
 *
 */
public class DNAFastaSymbolTableTest extends TestCase {

  /**
   */
  public DNAFastaSymbolTableTest(final String name) {
    super(name);
  }

  /**
   * Test method for {@link DNAFastaSymbolTable#scanResidue(int)}.
   */
  public final void testScanResidue() {
    DNAFastaSymbolTable t = new DNAFastaSymbolTable();
    assertEquals(DNA.A, t.scanResidue('A'));
    assertEquals(DNA.C, t.scanResidue('C'));
    assertEquals(DNA.G, t.scanResidue('G'));
    assertEquals(DNA.N, t.scanResidue('N'));
    assertEquals(DNA.T, t.scanResidue('T'));

    assertNull(t.scanResidue('J'));
    assertNull(t.scanResidue('O'));
  }

  /**
   * Test method for {@link DNAFastaSymbolTable#unknownResidue()}.
   */
  public final void testUnknownResidue() {
    DNAFastaSymbolTable table = new DNAFastaSymbolTable();
    assertEquals(DNA.N, table.unknownResidue());
  }

  /**
   * Test method for {@link DNAFastaSymbolTable#getSequenceType()}.
   */
  public final void testGetSequenceType() {
    DNAFastaSymbolTable t = new DNAFastaSymbolTable();
    assertEquals(SequenceType.DNA, t.getSequenceType());
  }

  /**
   * Test method for {@link DNAFastaSymbolTable#getAsciiToOrdinalTable()}
   */
  public final void testGetAsciiToOrdinalTable() {
    final DNAFastaSymbolTable t = new DNAFastaSymbolTable();
    final byte[] table = t.getAsciiToOrdinalTable();
    assertEquals(255, table.length);
    for (int i = 0; i < 255; i++) {
      final Residue r = t.scanResidue(i);
      if (r == null) {
        assertEquals((byte) 255, table[i]);
      } else {
        assertEquals((byte) r.ordinal(), table[i]);
      }
    }
  }

  /**
   * Test method for {@link DNAFastaSymbolTable#getOrdinalToAsciiTable()}
   */
  public final void testGetOrdinalToAsciiTable() {
    final DNAFastaSymbolTable t = new DNAFastaSymbolTable();
    final byte[] table = t.getOrdinalToAsciiTable();
    assertEquals(5, table.length);
    for (int k = 0; k < table.length; k++) {
      assertEquals((byte) DNA.values()[k].toString().charAt(0), table[k]);
    }
  }

}
