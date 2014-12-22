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
 * DNA implementation for FASTA symbol table.
 */
public class DNAFastaSymbolTable extends FastaSymbolTable {

  @Override
  public final byte[] getAsciiToOrdinalTable() {
    return getAsciiToOrdinalTable(this);
  }

  @Override
  public final byte[] getOrdinalToAsciiTable() {
    return getOrdinalToAsciiTable(DNA.values());
  }

  @Override
  public final Residue scanResidue(int code) {
    switch (code & (~32)) {
    case (int) 'A':
      return DNA.A;
    case (int) 'C':
      return DNA.C;
    case (int) 'G':
      return DNA.G;
    case (int) 'T':
    case (int) 'U':
      return DNA.T;
    case (int) 'N':
      return DNA.N;
    default:
    }
    return null;
  }

  @Override
  public final Residue unknownResidue() {
    return DNA.N;
  }

  @Override
  public final SequenceType getSequenceType() {
    return SequenceType.DNA;
  }
}

