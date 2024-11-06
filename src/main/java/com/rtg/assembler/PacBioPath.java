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
