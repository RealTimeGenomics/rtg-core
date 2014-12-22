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

/**
 */
public class AlignmentSection {
  long mContig;
  int mStartPosition;
  int mEndPosition;

  /**
   * record that a read aligns along <code>contig</code> from <code>startPosition</code> to <code>endPosition</code>
   * @param contig the contig id that this section aligns with
   * @param startPosition the start position within the contig (inclusive)
   * @param endPosition the end position within the contig (inclusive)
   */
  public AlignmentSection(long contig, int startPosition, int endPosition) {
    mContig = contig;
    mStartPosition = startPosition;
    mEndPosition = endPosition;
  }

  @Override
  public String toString() {
    return "AlignmentSection{" + "mContig=" + mContig + ", mStartPosition=" + mStartPosition + ", mEndPosition=" + mEndPosition + '}';
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }

    final AlignmentSection that = (AlignmentSection) o;

    return mContig == that.mContig
        && mEndPosition == that.mEndPosition
        && mStartPosition == that.mStartPosition;

  }

  @Override
  public int hashCode() {
    int result = (int) (mContig ^ (mContig >>> 32));
    result = 31 * result + mStartPosition;
    result = 31 * result + mEndPosition;
    return result;
  }
}
