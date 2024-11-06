/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package com.rtg.assembler;

import com.rtg.assembler.graph.Contig;
import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.mode.DNA;

/**
 * String representation of a <code>k-mer</code>.
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
      for (int i = start; i < end; ++i) {
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
    return new Factory();
  }

}
