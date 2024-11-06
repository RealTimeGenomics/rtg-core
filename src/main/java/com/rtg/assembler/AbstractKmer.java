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
    for (int i = 0; i < Math.min(this.length(), that.length()); ++i) {
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
    for (int i = 0; i < length(); ++i) {
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
    for (int i = 0; i < this.length(); ++i) {
      if (this.nt(i) != that.nt(i)) {
        return false;
      }
    }
    return true;
  }

  @Override
  public boolean integrity() {
    for (int i = 0; i < length(); ++i) {
      Exam.assertTrue(0 <= nt(i) && nt(i) <= 4);
    }
    return true;
  }

  @Override
  public void toString(StringBuilder sb) {
    //sb.append("k[" + length() + "]");
    final char[] valueChars = DNA.valueChars();
    for (int i = 0; i < length(); ++i) {
      final byte nt = nt(i);
      sb.append(valueChars[nt]);
    }
  }

}
