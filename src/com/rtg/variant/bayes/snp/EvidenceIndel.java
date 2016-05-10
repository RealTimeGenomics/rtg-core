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

package com.rtg.variant.bayes.snp;

import com.rtg.variant.bayes.EvidenceInterface;

/**
 */
public final class EvidenceIndel implements EvidenceInterface {

  /** Constant for a soft clip. */
  public static final int SOFT_CLIP_RIGHT = 3;

  /** Constant for a soft clip. */
  public static final int SOFT_CLIP_LEFT = 2;

  /** Constant for a deletion. */
  public static final int DELETE = 1;

  /** Constant for an insert. */
  public static final int INSERT = 0;

  /** Number of valid codes for <code>read</code> */
  public static final int DISTINCT_READS = 4;

  private final double mMapError;
  private final int mRead;

  private final int mMaxIndelLength;

  /**
   * @param mapError probability that the read doesn't map to this position.
   * @param indel 0 = insert, 1 = delete, 2 = soft clip left, 3 = soft clip right
   * @param maxIndelLength maximum length of an insert or delete as specified by cigar
   */
  public EvidenceIndel(final double mapError, int indel, int maxIndelLength) {
    mMapError = mapError;
    mRead = indel;
    mMaxIndelLength = maxIndelLength;
    assert mRead == INSERT || mRead == DELETE || mRead == SOFT_CLIP_LEFT || mRead == SOFT_CLIP_RIGHT;
  }

  @Override
  public double probability(int index) {
    throw new UnsupportedOperationException();
  }

  @Override
  public double pe() {
    throw new UnsupportedOperationException();
  }

  @Override
  public double error() {
    throw new UnsupportedOperationException();
  }

  @Override
  public double mapError() {
    return mMapError;
  }

  @Override
  public int read() {
    return mRead;
  }

  @Override
  public boolean isForward() {
    throw new UnsupportedOperationException();
  }

  //don't think these matter, since EvidenceIndel is merely used as a trigger...
  @Override
  public int getReadBasesLeft() {
    throw new UnsupportedOperationException();
  }

  @Override
  public int getReadBasesRight() {
    throw new UnsupportedOperationException();
  }

  @Override
  public boolean isReadPaired() {
    throw new UnsupportedOperationException();
  }

  @Override
  public boolean isMated() {
    throw new UnsupportedOperationException();
  }

  @Override
  public boolean isUnmapped() {
    throw new UnsupportedOperationException();
  }

  @Override
  public void setReadBasesLeft(int readBasesLeft) {
    throw new UnsupportedOperationException();
  }

  @Override
  public void setReadBasesRight(int readBaseRight) {
    throw new UnsupportedOperationException();
  }

  /**
   * @return maximum length of an insert or delete as specified by cigar
   */
  public int maxIndelLength() {
    return mMaxIndelLength;
  }

}
