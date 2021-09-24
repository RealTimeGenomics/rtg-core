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

import java.util.LinkedList;
import java.util.List;

/**
*/
class PacBioPath {
  int mReadStart;
  PacBioPath mPrevious;
  PartialAlignment mAlignment;
  boolean mIsPrefix = false;
  final int mScore;
  final long mFirstContig;
  final int mContigStart;
  boolean mIsDupe;
  PacBioPath(PacBioPath prev, PartialAlignment alignment) {
    if (prev == null) {
      mReadStart = alignment.getReadStart();
      mScore = alignment.getAlignmentScore();
      mFirstContig = alignment.getContig();
      mContigStart = alignment.getContigStart();
      mIsDupe = false;
    } else {
      mReadStart = prev.mReadStart;
      mScore = prev.mScore + alignment.getAlignmentScore();
      mContigStart = prev.mContigStart;
      mFirstContig = prev.mFirstContig;
      mIsDupe = prev.mIsDupe;
    }
    mPrevious = prev;
    mAlignment = alignment;
  }

  int score() {
    return mScore;
  }
  List<Long> toPath() {
    PacBioPath current = this;
    final List<Long> path = new LinkedList<>();
    while (current != null) {
      path.add(0, current.mAlignment.getContig());
      current = current.mPrevious;
    }
    return path;

  }

  @Override
  public String toString() {
    final StringBuilder pathString = new StringBuilder();
    int score = 0;
    int finalReadStart = 0;
    PacBioPath current = this;
    while (current != null) {
      pathString.insert(0, current.mAlignment.getContig());
      pathString.insert(0, ", ");
      score += current.mAlignment.getAlignmentScore();
      finalReadStart = current.mAlignment.getReadStart();
      current = current.mPrevious;
    }
    return pathString + " isDupe=" + mIsDupe + " score=" + score + " readLength=" + (this.mAlignment.getReadEnd() - finalReadStart) + "[" + mReadStart + "-" + this.mAlignment.getReadEnd() + "]";
  }
}
