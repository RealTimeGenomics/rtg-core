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

package com.rtg.assembler.graph.implementation;

import com.rtg.assembler.graph.Contig;
import com.rtg.mode.DNA;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.mode.DNARange;

/**
 */
public class ContigString extends IntegralAbstract implements Contig {

  private final byte[] mNt;

  /**
   * @param nt string containing the nucleotides in <code>NACGT</code> format.
   */
  public ContigString(final String nt) {
    mNt = new byte[nt.length()];
    for (int i = 0; i < nt.length(); ++i) {
      mNt[i] = DNARange.RANGE.valueOf(nt.charAt(i));
    }
  }

  @Override
  public boolean integrity() {
    return true;
  }

  @Override
  public int length() {
    return mNt.length;
  }

  @Override
  public byte nt(int index) {
    return mNt[index];
  }

  /**
   * @param c the contig to convert to a string
   * @return a string representation of the bases within this contig
   */
  public static String contigSequenceString(Contig c) {
    return contigSubString(c, 0, c.length());
  }

  /**
   * @param c the contig to convert to a string
   * @param start first base position
   * @param end exclusive base position
   * @return a string representation of the bases within this contig
   */
  public static String contigSubString(Contig c, int start, int end) {
    final StringBuilder sb = new StringBuilder();
    for (int i = start; i < end; ++i) {
      sb.append(DNA.valueChars()[c.nt(i)]);
    }
    return sb.toString();
  }
}
