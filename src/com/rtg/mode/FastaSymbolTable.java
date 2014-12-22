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

import com.reeltwo.jumble.annotations.TestClass;


/**
 * Interface for FASTA symbol table.
 *
 */
@TestClass("com.rtg.mode.DNAFastaSymbolTableTest")
public abstract class FastaSymbolTable {
  /**
   * Takes an <code>int</code> code and returns the associated residue
   * @param code byte read from FASTA file
   * @return residue representing code
   */
  public abstract Residue scanResidue(final int code);

  /**
   * Returns the residue to use for unknown values
   * @return the residue
   */
  public abstract Residue unknownResidue();

  /**
   * Returns the sequence type for current table
   * @return the sequence type
   */
  public abstract SequenceType getSequenceType();

  /**
   * Returns a table that can be used to map from the internal code
   * (ordinal value) to an ASCII character representation. The table
   * will have one entry per possible residue.
   * @return the mapping table.
   */
  public abstract byte[] getOrdinalToAsciiTable();

  /**
   * Returns a table that can be used to map from the character
   * representation to an internal code (the ordinal of the
   * residue). Any unrecognized characters will map in this table to
   * value 255.
   * @return the mapping table.
   */
  public abstract byte[] getAsciiToOrdinalTable();

  // Creates map from ascii representation to internal byte code
  static byte[] getAsciiToOrdinalTable(FastaSymbolTable table) {
    final byte[] mapping = new byte[255];
    Residue r;
    for (int i = 0; i < 255; i++) {
      r = table.scanResidue(i);
      if (r == null) {
        mapping[i] = (byte) 255;
      } else {
        mapping[i] = (byte) r.ordinal();
      }
    }
    return mapping;
  }

  // Creates map from internal byte code to ascii representation
  static byte[] getOrdinalToAsciiTable(Residue[] residues) {
    final byte[] mapping = new byte[residues.length];
    for (int k = 0; k < mapping.length; k++) {
      mapping[k] = (byte) residues[k].toString().charAt(0);
    }
    return mapping;
  }
}
