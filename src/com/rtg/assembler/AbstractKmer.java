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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.mode.DNA;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Implement equals and hash for all Kmers.
 */
@TestClass("com.rtg.assembler.StringKmerTest")
public abstract class AbstractKmer extends IntegralAbstract implements Kmer {

  @Override
  public int compareTo(Kmer that) {
    for (int i = 0; i < Math.min(this.length(), that.length()); i++) {
      final byte thisNt = this.nt(i);
      final byte thatNt = that.nt(i);
      final int delta = thisNt - thatNt;
      if (delta != 0) {
        return delta;
      }
    }
    return this.length() - that.length();
  }

  @Override
  public int hashCode() {
    int hash = 31;
    for (int i = 0; i < length(); i++) {
      hash = hash * 1000000007 + nt(i);
    }
    return hash;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || !(o instanceof Kmer)) {
      return false;
    }

    final Kmer that = (Kmer) o;
    if (this.length() != that.length()) {
      return false;
    }
    for (int i = 0; i < this.length(); i++) {
      if (this.nt(i) != that.nt(i)) {
        return false;
      }
    }
    return true;
  }

  @Override
  public boolean integrity() {
    for (int i = 0; i < length(); i++) {
      Exam.assertTrue(0 <= nt(i) && nt(i) <= 4);
    }
    return true;
  }

  @Override
  public void toString(StringBuilder sb) {
    //sb.append("k[" + length() + "]");
    final char[] valueChars = DNA.valueChars();
    for (int i = 0; i < length(); i++) {
      final byte nt = nt(i);
      sb.append(valueChars[nt]);
    }
  }

}
