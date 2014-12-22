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

import java.util.HashMap;
import java.util.Iterator;

import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;

/**
 *         Date: 12/04/12
 *         Time: 5:38 PM
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
    for (int i = 0; i < bases.length(); i++) {
      b[i] = bases.baseAt(i);
    }
    return b;
  }
  private static String sequenceBasesToString(SequenceBases bases) {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < bases.length(); i++) {
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

    public Mutations(SequenceBases original, boolean all, int start, int end) {
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
      mIndex++;
      while (mIndex == mMutant.baseAt(mPosition)) {
        mIndex++;
      }
      if (mIndex > T) {
        mIndex = A;
        mPosition++;
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
