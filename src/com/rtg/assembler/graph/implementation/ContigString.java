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

package com.rtg.assembler.graph.implementation;

import com.rtg.assembler.graph.Contig;
import com.rtg.mode.DNA;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.dna.DNARange;

/**
 */
public class ContigString extends IntegralAbstract implements Contig {

  private final byte[] mNt;

  /**
   * @param nt string containing the nucleotides in <code>NACGT</code> format.
   */
  public ContigString(final String nt) {
    mNt = new byte[nt.length()];
    for (int i = 0; i < nt.length(); i++) {
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
    for (int i = start; i < end; i++) {
      sb.append(DNA.valueChars()[c.nt(i)]);
    }
    return sb.toString();
  }
}
