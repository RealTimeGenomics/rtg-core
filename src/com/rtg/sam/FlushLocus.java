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

package com.rtg.sam;

import com.rtg.util.Utils;

/**
 */
class FlushLocus {
  int mStart;
  int mEnd;

  public FlushLocus(int start, int end) {
    mStart = start;
    mEnd = end;
  }

  public boolean isJoinable(FlushLocus other) {
    return other.mStart <= mEnd && other.mEnd >= mStart;
  }

  public void join(FlushLocus other) {
    mStart = Math.min(mStart, other.mStart);
    mEnd = Math.max(mEnd, other.mEnd);
  }

  @Override
  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    if (obj == this) {
      return true;
    }
    final FlushLocus that = (FlushLocus) obj;
    return this.mStart == that.mStart && this.mEnd == that.mEnd;
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(mStart, mEnd);
  }

  @Override
  public String toString() {
    return "[" + mStart + "," + mEnd + ")";
  }
}
