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

import java.util.HashMap;
import java.util.Iterator;

import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;

/**
 */
public class CorrectingMutator {
  static final byte A = (byte) DNA.A.ordinal();
  static final byte T = (byte) DNA.T.ordinal();


  CorrectingMutator() {
  }
  Iterable<SequenceBases> getMutations(SequenceBases original, int start, int end) {
    return new Mutations(original, true, start, end);
  }

  /**
   * Interface providing an override-able way of retrieving the bases in a sequence. Useful for representing error
   * corrected reads without modifying an underlying array.
   */
  public interface SequenceBases {
    /**
     * @param position within the sequence
     * @return the base at <code>position</code> in the sequence
     */
    byte baseAt(int position);
    /**
     * @return total sequence length
     */
    int length();
  }

  /**
   * Represents a read without mutations. Simply passes the call through to an array lookup
   */
  public static class BaseRead implements SequenceBases {
    final byte[] mOriginal;

    /**
     * @param read the read read bases
     */
    public BaseRead(byte[] read) {
      mOriginal = read;
    }

    @Override
    public byte baseAt(int position) {
      return mOriginal[position];
    }

    @Override
    public int length() {
      return mOriginal.length;
    }

    @Override
    public String toString() {
      return sequenceBasesToString(this);
    }
  }
  static byte[] sequenceBasesToBytes(SequenceBases bases) {
    final byte[] b = new byte[bases.length()];
    for (int i = 0; i < bases.length(); ++i) {
      b[i] = bases.baseAt(i);
    }
    return b;
  }
  private static String sequenceBasesToString(SequenceBases bases) {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < bases.length(); ++i) {
      sb.append(DnaUtils.getBase(bases.baseAt(i)));
    }
    return sb.toString();
  }

  /**
   * Represents a read with single some nucleotide errors corrected
   */
  public static class CorrectedRead implements SequenceBases {
    HashMap<Integer, Byte> mMutations = new HashMap<>();
    SequenceBases mOriginal;

    /**
     * @param original the original read bases
     */
    public CorrectedRead(SequenceBases original) {
      mOriginal = original;
    }

    @Override
    public byte baseAt(int position) {
      if (mMutations.containsKey(position)) {
        return mMutations.get(position);
      }
      return mOriginal.baseAt(position);
    }

    /**
     * Add a nucleotide correction
     * @param position position to correct
     * @param base the corrected base
     */
    public void correct(int position, byte base) {
      mMutations.put(position, base);
    }

    @Override
    public int length() {
      return mOriginal.length();
    }
    @Override
    public String toString() {
      return sequenceBasesToString(this);
    }
  }

  private static class Mutations implements Iterable<SequenceBases> {
    SequenceBases mOriginal;
    boolean mAll;
    int mStart;
    int mEnd;

    Mutations(SequenceBases original, boolean all, int start, int end) {
      mOriginal = original;
      mAll = all;
      mStart = start;
      mEnd = end;
    }

    @Override
    public Iterator<SequenceBases> iterator() {
      if (mAll) {
        return new AllMutations(mOriginal, mStart, mEnd);
      } else {
        throw new UnsupportedOperationException();
      }
    }
  }

  private static class AllMutations implements Iterator<SequenceBases> {
    SequenceBases mMutant;
    int mPosition;
    byte mIndex;
    int mEnd;
    AllMutations(SequenceBases original, int start, int end) {
      mMutant = original;
      mPosition = start;
      mIndex = (byte) DNA.N.ordinal();
      mEnd = end;
    }

    @Override
    public boolean hasNext() {
      final byte lastBase = mMutant.baseAt(mEnd - 1);
      return mPosition < mEnd - 1 || mIndex < T && lastBase < T || lastBase == T && mIndex < T - 1;
    }

    @Override
    public SequenceBases next() {
      ++mIndex;
      while (mIndex == mMutant.baseAt(mPosition)) {
        ++mIndex;
      }
      if (mIndex > T) {
        mIndex = A;
        ++mPosition;
      }
      final CorrectedRead corrected = new CorrectedRead(mMutant);
      corrected.correct(mPosition, mIndex);
      return corrected;
    }

    @Override
    public void remove() {
      throw new UnsupportedOperationException();
    }
  }
}
