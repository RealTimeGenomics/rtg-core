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
    if (this.mEnd < that.mEnd) {
      return -1;
    }
    if (this.mEnd > that.mEnd) {
      return +1;
    }
    return 0;
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
