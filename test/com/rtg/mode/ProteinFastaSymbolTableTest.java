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
public class ProteinFastaSymbolTableTest extends TestCase {

  /**
   */
  public ProteinFastaSymbolTableTest(final String name) {
    super(name);
  }

  /**
   * Test method for {@link ProteinFastaSymbolTable#scanResidue(int)}.
   */
  public final void testScanResidue() {
    ProteinFastaSymbolTable t = new ProteinFastaSymbolTable();
    assertEquals(Protein.A, t.scanResidue('A'));
    assertEquals(Protein.C, t.scanResidue('C'));
    assertEquals(Protein.D, t.scanResidue('D'));
    assertEquals(Protein.E, t.scanResidue('E'));
    assertEquals(Protein.F, t.scanResidue('F'));
    assertEquals(Protein.G, t.scanResidue('G'));
    assertEquals(Protein.H, t.scanResidue('H'));
    assertEquals(Protein.I, t.scanResidue('I'));
    assertEquals(Protein.K, t.scanResidue('K'));
    assertEquals(Protein.L, t.scanResidue('L'));
    assertEquals(Protein.M, t.scanResidue('M'));
    assertEquals(Protein.N, t.scanResidue('N'));
    assertEquals(Protein.P, t.scanResidue('P'));
    assertEquals(Protein.Q, t.scanResidue('Q'));
    assertEquals(Protein.R, t.scanResidue('R'));
    assertEquals(Protein.S, t.scanResidue('S'));
    assertEquals(Protein.T, t.scanResidue('T'));
    assertEquals(Protein.V, t.scanResidue('V'));
    assertEquals(Protein.W, t.scanResidue('W'));
    assertEquals(Protein.Y, t.scanResidue('Y'));
    assertEquals(Protein.X, t.scanResidue('X'));
    assertEquals(Protein.X, t.scanResidue('Z'));
    assertEquals(Protein.X, t.scanResidue('B'));

    assertNull(t.scanResidue('J'));
    assertNull(t.scanResidue('O'));
    assertNull(t.scanResidue('U'));

    assertEquals(Protein.STOP, t.scanResidue('*'));
  }

  /**
   * Test method for {@link ProteinFastaSymbolTable#unknownResidue()}.
   */
  public final void testUnknownResidue() {
    ProteinFastaSymbolTable table = new ProteinFastaSymbolTable();
    assertEquals(Protein.X, table.unknownResidue());
  }

  /**
   * Test method for {@link ProteinFastaSymbolTable#getSequenceType()}.
   */
  public final void testGetSequenceType() {
    ProteinFastaSymbolTable t = new ProteinFastaSymbolTable();
    assertEquals(SequenceType.PROTEIN, t.getSequenceType());
  }

  /**
   * Test method for {@link ProteinFastaSymbolTable#getAsciiToOrdinalTable()}
   */
  public final void testGetAsciiToOrdinalTable() {
    final ProteinFastaSymbolTable t = new ProteinFastaSymbolTable();
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
   * Test method for {@link ProteinFastaSymbolTable#getOrdinalToAsciiTable()}
   */
  public final void testGetOrdinalToAsciiTable() {
    final ProteinFastaSymbolTable t = new ProteinFastaSymbolTable();
    final byte[] table = t.getOrdinalToAsciiTable();
    assertEquals(22, table.length);
    for (int k = 0; k < table.length; k++) {
      assertEquals((byte) Protein.values()[k].toString().charAt(0), table[k]);
    }
  }

}
