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
package com.rtg.position.output;

import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Holds a start and end position corresponding to an exact segment match.
 */
public class Segment extends IntegralAbstract implements Comparable<Segment> {

  private static final int NULL = -1;

  private int mStart = NULL;

  private int mEnd = NULL;

  private int mSeqId = NULL;

  /**
   * Create initial valid state.
   * @param seqid identifier of the sequence.
   * @param start position of first residue in segment (counting from 0).
   * @param end position of last residue in segment (counting from 0).
   */
  public void initialize(final int seqid, final int start, final int end) {
    assert start >= 0 && end >= start; // : "start=" + start + "end=" + end;
    assert mStart == NULL && mEnd == NULL && mSeqId == NULL;
    mSeqId = seqid;
    mStart = start;
    mEnd = end;
    integrity();
  }

  /**
   * Clear the state so initialize can be called later (same as initial state).
   */
  public void clear() {
    mSeqId = NULL;
    mStart = NULL;
    mEnd = NULL;
  }

  /**
   * Extend segment given a new end position.
   * Must be later than the current end position.
   * @param end last residue position for the match (counts from 0).
   */
  public void extend(final int end) {
    assert end > mEnd;
    mEnd = end;
  }

  /**
   * Test if the segment is empty.
   * @return true iff the segment is empty.
   */
  public boolean isEmpty() {
    return mStart == NULL && mEnd == NULL && mSeqId == NULL;
  }

  /**
   * Get the start position.
   * @return the start position.
   */
  public int start() {
    return mStart;
  }

  /**
   * Get the end position.
   * @return the end position. 0 based inclusive
   */
  public int end() {
    return mEnd;
  }

  /**
   * Get the sequence id.
   * @return the sequence id.
   */
  public int seqId() {
    return mSeqId;
  }

  @Override
  public boolean equals(final Object obj) {
    if (this == obj) {
      return true;
    }
    if (!(obj instanceof Segment)) {
      return false;
    }
    final Segment that = (Segment) obj;
    return this.mSeqId == that.mSeqId && this.mStart == that.mStart && this.mEnd == that.mEnd;
  }


  @Override
  public int hashCode() {
    final long x = Utils.pairHash(mSeqId, Utils.pairHash(mStart, mEnd));
    return (int) x;
  }

  @Override
  public int compareTo(final Segment that) {
    if (this == that) {
      return 0;
    }
    if (this.mSeqId < that.mSeqId) {
      return -1;
    }
    if (this.mSeqId > that.mSeqId) {
      return +1;
    }
    if (this.mStart < that.mStart) {
      return -1;
    }
    if (this.mStart > that.mStart) {
      return +1;
    }
    return Integer.compare(this.mEnd, that.mEnd);
  }

  @Override
  public void toString(final StringBuilder sb) {
    if (mStart == NULL) {
      sb.append("empty");
    } else {
      sb.append(mSeqId).append(":").append(mStart).append("..").append(mEnd);
    }
  }

  @Override
  public boolean integrity() {
    if (mStart == NULL) {
      Exam.assertEquals(NULL, mEnd);
      Exam.assertEquals(NULL, mSeqId);
    } else {
      Exam.assertTrue(mSeqId >= 0);
      Exam.assertTrue(mStart >= 0);
      Exam.assertTrue(mEnd >= mStart);
    }
    return true;
  }

}
