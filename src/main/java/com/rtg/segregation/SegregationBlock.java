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

package com.rtg.segregation;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Represents a consecutive sequence of segregated calls that all have compatible phasing with each other.
 */
public class SegregationBlock extends IntegralAbstract {

  private final String mSeq;
  private final int mStart;
  private int mEnd;
  private int mCount = 1;
  private PatternArray mPatterns;
  private final boolean mIsXlike;

  /**
   * @param family family info
   */
  public SegregationBlock(final FamilyGt family) {
    mSeq = family.seq();
    mStart = family.posn();
    mEnd = family.posn();
    mPatterns = family.pattern();
    mIsXlike = family.isXLike();
    assert integrity();
  }

  /**
   * @return the pattern array.
   */
  PatternArray patterns() {
    return mPatterns;
  }

  /**
   * Extend the block if it can be done and keep consistent phasing. Update the block phasing pattern if successful.
   * @param family containing a set of calls.
   * @return true iff the block was successfully extended.
   */
  boolean extend(final FamilyGt family) {
    final int length = family.length();
    assert length == mPatterns.length();
    assert family.posn() > mEnd;
    assert mSeq.equals(family.seq());
    final PatternArray pattern = family.pattern();
    for (int f = 0; f < Pattern.NUMBER_FLIPS; ++f) {
      final PatternArray pat = mPatterns.flipIntersect(pattern, f);
      if (pat != null) {
        mPatterns = pat;
        mEnd = family.posn();
        ++mCount;
        return true;
      }
    }
    return false;
  }

  int count() {
    return mCount;
  }

  int start() {
    return mStart;
  }

  int end() {
    return mEnd;
  }

  boolean isXLike() {
    return mIsXlike;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append(mSeq);
    sb.append(" ").append(mStart).append(" - ").append(mEnd);
    sb.append(" ").append(mCount);
    sb.append(" fa ").append(patterns().faString());
    sb.append(" mo ").append(patterns().moString());
  }

  @Override
  public final boolean integrity() {
    Exam.assertTrue(1 <= mStart);
    Exam.assertTrue(mStart <= mEnd);
    Exam.assertNotNull(mPatterns);
    Exam.assertTrue(mPatterns.length() > 1);
    Exam.assertTrue(mCount >= 1);
    if (mEnd > mStart) {
      Exam.assertTrue(mCount > 1);
    }
    return true;
  }
}
