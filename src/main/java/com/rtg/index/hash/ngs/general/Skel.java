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
package com.rtg.index.hash.ngs.general;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Holds information about the contiguous bits in a mask.
 */
public final class Skel extends IntegralAbstract {

  private static final long LONG_BITS = 64;

  /** The bit position of the highest order bit (counting from 0). */
  private final int mPosition;

  /** The number of bits. */
  private final int mLength;

  /** The number of bits to right shift to get to all bits in the low order of the long. */
  private final int mRtShift;

  /** The number of bits to left shift after masking to get into correct position. */
  private final int mLeftShift;


  /**
   * @param position the position
   * @param length the length
   * @param finalPosition the final position
   */
  public Skel(final int position, final int length, final int finalPosition) {
    //System.err.println("Skel position=" + position + " length=" + length + " finalPosition=" + finalPosition);
    mLength = length;
    mPosition = position;
    mRtShift = mPosition - mLength + 1;
    mLeftShift = finalPosition - mLength + 1;
    integrity();
  }

  /**
   * Extract the bits.
   * @param v original value the bits are to be extracted from.
   * @return the final extracted bits.
   */
  public long mask(final long v) {
    final long t = v >>> mRtShift;
    final long u = t & (mLength == 64 ? -1L : ((1L << mLength) - 1L));
    return u << mLeftShift;
  }

  /**
   * Extract the bits making allowance for an indel tweak.
   * @param v original value the bits are to be extracted from.
   * @param tweak the increment to the start position used to allow for for indels.
   * @return the final extracted bits.
   */
  public long mask(final long v, final int tweak) {
    assert mRtShift + tweak >= 0;
    final long t = v >>> (mRtShift + tweak);
    final long u = t & (mLength == 64 ? -1L : ((1L << mLength) - 1L));
    return u << mLeftShift;
  }

  /**
   * Get the number of adjacent positions covered.
   * @return the number of adjacent positions covered.
   */
  int length() {
    return mLength;
  }

  /**
   * Get the highest position covered (based from 0).
   * @return the highest position covered.
   */
  int position() {
    return mPosition;
  }

  /**
   * Get the highest position after masking (based from 0).
   * @return the highest position after masking.
   */
  int finalPosition() {
    return mLeftShift + mLength - 1;
  }

  @Override
  public void toString(final StringBuilder sb) {
    //sb.append("0b").append(Utils.toBits(((1L << mLength) - 1) << mRtShift, Long.SIZE)).append(' '); // More visual representation of the skeleton blocks
    sb.append("[").append(mPosition).append("...").append(mPosition - mLength + 1).append("]<<").append(mLeftShift);
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(0 <= mPosition && mPosition < LONG_BITS);
    Exam.assertTrue(mPosition >= (mLength - 1));
    Exam.assertTrue(mLength >= 1);
    Exam.assertTrue(0 <= mRtShift && mRtShift < LONG_BITS);
    Exam.assertTrue(this.toString() + " leftShift=" + mLeftShift + " rtShift=" + mRtShift, 0 <= mLeftShift && mLeftShift <= mRtShift);
    return true;
  }

}
