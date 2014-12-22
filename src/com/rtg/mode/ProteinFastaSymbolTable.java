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


/**
 * Protein implementation for FASTA symbol table.
 */
public class ProteinFastaSymbolTable extends FastaSymbolTable {

  @Override
  public final byte[] getAsciiToOrdinalTable() {
    return getAsciiToOrdinalTable(this);
  }

  @Override
  public final byte[] getOrdinalToAsciiTable() {
    return getOrdinalToAsciiTable(Protein.values());
  }


  @Override
  public Residue scanResidue(final int code) {
    // This is done first, so that ~32 does not interfere!
    if (code == '*') {
      return Protein.STOP;
    }
    switch (code & (~32)) {
    case (int) 'A':
      return Protein.A;
    case (int) 'C':
      return Protein.C;
    case (int) 'D':
      return Protein.D;
    case (int) 'E':
      return Protein.E;
    case (int) 'F':
      return Protein.F;
    case (int) 'G':
      return Protein.G;
    case (int) 'H':
      return Protein.H;
    case (int) 'I':
      return Protein.I;
    case (int) 'K':
      return Protein.K;
    case (int) 'L':
      return Protein.L;
    case (int) 'M':
      return Protein.M;
    case (int) 'N':
      return Protein.N;
    case (int) 'P':
      return Protein.P;
    case (int) 'Q':
      return Protein.Q;
    case (int) 'R':
      return Protein.R;
    case (int) 'S':
      return Protein.S;
    case (int) 'T':
      return Protein.T;
    case (int) 'V':
      return Protein.V;
    case (int) 'W':
      return Protein.W;
    case (int) 'Y':
      return Protein.Y;
    case (int) 'X': //unknown residue
    case (int) 'Z': //Z has 2 options
    case (int) 'B': //B has 2 options
      return Protein.X;
    default:
      return null;
    }
  }
  @Override
  public Residue unknownResidue() {
    return Protein.X;
  }
  @Override
  public SequenceType getSequenceType() {
    return SequenceType.PROTEIN;
  }
}
