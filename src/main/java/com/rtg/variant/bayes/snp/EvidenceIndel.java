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

package com.rtg.variant.bayes.snp;

import com.rtg.variant.bayes.EvidenceInterface;

/**
 * Evidence for an indel event.
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

  private final int mOperationLength;

  /**
   * @param mapError probability that the read doesn't map to this position.
   * @param indel 0 = insert, 1 = delete, 2 = soft clip left, 3 = soft clip right
   * @param operationLength maximum length of an insert or delete as specified by cigar
   */
  public EvidenceIndel(final double mapError, int indel, int operationLength) {
    mMapError = mapError;
    mRead = indel;
    mOperationLength = operationLength;
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
  public boolean isFirst() {
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
  public int operationLength() {
    return mOperationLength;
  }

}
