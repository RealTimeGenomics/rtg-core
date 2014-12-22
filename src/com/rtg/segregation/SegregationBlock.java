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
    for (int f = 0; f < Pattern.NUMBER_FLIPS; f++) {
      final PatternArray pat = mPatterns.flipIntersect(pattern, f);
      if (pat != null) {
        mPatterns = pat;
        mEnd = family.posn();
        mCount++;
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
  public boolean integrity() {
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
