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

import com.rtg.util.CompareHelper;
import com.rtg.util.Utils;

/**
*/
class PartialAlignment implements Comparable<PartialAlignment> {
  private final int mAlignmentScore;
  private final int mReadStart;
  private final int mReadEnd;
  private final int mContigStart;
  private final int mContigEnd;
  private final long mContig;

  PartialAlignment(int alignmentScore, int readStart, int readEnd, long contig, int contigStart, int contigEnd) {
    mAlignmentScore = alignmentScore;
    mReadStart = readStart;
    mReadEnd = readEnd;
    mContig = contig;
    mContigStart = contigStart;
    mContigEnd = contigEnd;
  }

  public long getContig() {
    return mContig;
  }

  public int getAlignmentScore() {
    return mAlignmentScore;
  }

  public int getReadStart() {
    return mReadStart;
  }

  public int getReadEnd() {
    return mReadEnd;
  }

  public int getContigStart() {
    return mContigStart;
  }

  public int getContigEnd() {
    return mContigEnd;
  }

  @Override
  public int compareTo(PartialAlignment other) {
    return new CompareHelper()
        .compare(getReadStart(), other.getReadStart())
        .compare(getReadEnd(), other.getReadEnd())
        .compare(getContig(), other.getContig())
        .compare(getContigStart(), other.getContigStart())
        .compare(getContigEnd(), other.getContigEnd())
        .compare(getAlignmentScore(), other.getAlignmentScore())
        .result();
  }

  @Override
  public String toString() {
    return "contig=" + getContig() + " [" + getContigStart() + "," + getContigEnd() + "] read [" + getReadStart() + "," + getReadEnd() + "] score=" + getAlignmentScore();
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }

    final PartialAlignment that = (PartialAlignment) o;
    return this.compareTo(that) == 0;
  }

  @Override
  public int hashCode() {
    return Utils.pairHashContinuous(getAlignmentScore(), getReadStart(), getReadEnd(), getContigStart(), getContigEnd(), Long.valueOf(getContig()).hashCode());
  }
}
