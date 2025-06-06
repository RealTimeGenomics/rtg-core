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
    return 31 * result + mEndPosition;
  }
}
