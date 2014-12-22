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
