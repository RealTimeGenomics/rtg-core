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

package com.rtg.variant.eval;

import java.util.Arrays;

import com.rtg.mode.DnaUtils;
import com.rtg.util.intervals.SequenceNameLocus;
import com.rtg.util.Utils;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class MockVariant extends IntegralAbstract implements Variant {

  private final int mStart;
  private final int mEnd;
  private final byte[] mPlus;
  private final byte[] mMinus;
  private final boolean mPhased;

  /**
   * @param start one-based start position of mutation
   * @param end one-based end position of mutation
   * @param plus nucleotides on the plus strand
   * @param minus nucleotides on the minus strand
   * @param phased does this call have phasing information
   */
  public MockVariant(int start, int end, byte[] plus, byte[] minus, boolean phased) {
    super();
    mStart = start - 1;
    mEnd = end - 1;
    mPlus = Arrays.copyOf(plus, plus.length);
    if (minus != null) {
      mMinus = Arrays.copyOf(minus, minus.length);
    } else {
      mMinus = null;
    }
    mPhased = phased;
  }
  /**
   * Assumes not phased
   * @param start one-based start position of mutation
   * @param end one-based end position of mutation
   * @param plus nucleotides on the plus strand
   * @param minus nucleotides on the minus strand
   */
  public MockVariant(int start, int end, byte[] plus, byte[] minus) {
    this(start, end, plus, minus, false);
  }

  @Override
  public int getEnd() {
    return mEnd;
  }

  @Override
  public boolean overlaps(SequenceNameLocus other) {
    throw new UnsupportedOperationException();
  }

  @Override
  public boolean contains(String sequence, int pos) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int getLength() {
    return getEnd() - getStart();
  }

  @Override
  public byte[] nt(boolean strand) {
    return strand ? ntAlleleA() : ntAlleleB();
  }

  @Override
  public byte[] ntAlleleB() {
    if (mMinus != null) {
      return Arrays.copyOf(mMinus, mMinus.length);
    } else {
      return null;
    }
  }

  @Override
  public byte[] ntAlleleA() {
    return Arrays.copyOf(mPlus, mPlus.length);
  }

  @Override
  public int getStart() {
    return mStart;
  }

  @Override
  public boolean integrity() {
    return true;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append(getStart() + 1).append(":").append(getEnd() + 1).append(" ");
    sb.append(DnaUtils.bytesToSequenceIncCG(mPlus));
    if (mMinus != null) {
      sb.append(":").append(DnaUtils.bytesToSequenceIncCG(mMinus));
    }
  }

  @Override
  public boolean isPhased() {
    return mPhased;
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(mStart, mEnd);
  }

  @Override
  public boolean equals(final Object o2) {
    if (this == o2) {
      return true;
    }
    if (o2 == null) {
      return false;
    }
    if (!(o2 instanceof MockVariant)) {
      return false;
    }
    final MockVariant other = (MockVariant) o2;
    return mStart == other.mStart && mEnd == other.mEnd;
  }

  @Override
  public String getSequenceName() {
    return "";
  }
}
