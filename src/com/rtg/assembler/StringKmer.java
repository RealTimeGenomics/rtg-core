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

package com.rtg.assembler;

import com.rtg.assembler.graph.Contig;
import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.mode.DNA;

/**
 *         Date: 27/04/12
 *         Time: 2:04 PM
 */
class StringKmer extends AbstractKmer {
  final String mCode;

  StringKmer(String code) {
    mCode = code;
  }

  @Override
  public Kmer predecessor(byte nt) {
    return new StringKmer(DNA.valueChars()[nt] + mCode.substring(0, length() - 1));
  }

  @Override
  public Kmer successor(byte nt) {
    return new StringKmer(mCode.substring(1) + DNA.valueChars()[nt]);
  }

  @Override
  public int length() {
    return mCode.length();
  }

  @Override
  public byte nt(int pos) {
    return (byte) DNA.getDNA(mCode.charAt(pos));
  }

  @Override
  public Kmer minimalKmer() {
    final String minimal = DeBruijnGraphBuilder.minimalVersion(mCode);
    if (minimal.compareTo(mCode) < 0) {
      return new StringKmer(minimal);
    } else {
      return this;
    }
  }
  @Override
  public Kmer reverse() {
    return new StringKmer(DeBruijnGraphBuilder.complement(mCode));
  }

  private static class Factory implements KmerFactory {
    @Override
    public Kmer make(byte[] kmer, int start, int end) {
      final StringBuilder sb = new StringBuilder();
      for (int i = start; i < end; i++) {
        sb.append(DNA.valueChars()[kmer[i]]);
      }
      return new StringKmer(sb.toString());
    }
    @Override
    public Kmer make(Contig contig, int start, int end) {
      return new StringKmer(ContigString.contigSubString(contig, start, end));
    }
  }

  static KmerFactory factory() {
    //Diagnostic.userLog("Using string kmer factory");
    return new Factory();
  }

}
