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
package com.rtg.index.queue;

import com.rtg.util.array.ExtensibleIndex;
import com.rtg.util.array.intindex.IntChunks;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;


/**
 */
public class QueueIterator extends IntegralAbstract {

  private final ExtensibleIndex mMemory;

  private final long mRadixShifted;

  private long mEnd;

  private long mCurrent;

  //last position in this block (exclusive)
  private long mLast;

  private long mNext;

  private long mHash;

  private long mId;

  private boolean mOK;

  QueueIterator(final ExtensibleIndex memory, final long start, final long radixShifted) {
    mMemory = memory;
    mRadixShifted = radixShifted;
    step(start);
  }

  private long getFree(final long pos) {
    if (mMemory instanceof IntChunks) {
      return (mMemory.getSigned(pos - 1) << 32) | mMemory.get(pos);
    } else {
      assert mMemory.getSigned(pos - 1) == 0;
      return mMemory.getSigned(pos);
    }
  }

  private void step(final long start) {
    mEnd = start;
    final long length = mMemory.getSigned(start);
    final long s1 = getFree(start - 1);
    if (length < 0) {
      mCurrent = mEnd + length + 1;
      mLast = s1;
      mNext = -1;
    } else {
      mCurrent = mEnd - length + 1;
      mLast = mEnd - 2;
      mNext = s1;
    }
  }

  /**
   * Step to the next hash, id pair.
   * @return true iff a new valid pair is available.
   */
  boolean next() {
    if (mCurrent == mLast) {
      if (mNext < 0) {
        mOK = false;
        return false;
      }
      step(mNext);
    }
    mHash = mMemory.get(mCurrent) | mRadixShifted;
    ++mCurrent;

    if (mCurrent == mLast) {
      assert mNext >= 0;
      step(mNext);
    }
    mId = mMemory.get(mCurrent);
    ++mCurrent;
    mOK = true;
    return true;
  }

  /**
   * Get the hash value.
   * @return  the hash value.
   */
  long hash() {
    if (!mOK) {
      throw new IllegalStateException();
    }
    return mHash;
  }

  /**
   * Get the id value.
   * @return  the id value.
   */
  long id() {
    if (!mOK) {
      throw new IllegalStateException();
    }
    return mId;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("QueueIterator:");
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mMemory);
    final long memLength = mMemory.length();
    Exam.assertTrue(0 <= mCurrent && mCurrent <= mLast);
    Exam.assertTrue(mLast <= mEnd && mEnd < memLength);
    Exam.assertTrue(mNext == -1 || (mEnd < mNext && mNext < memLength));
    return true;
  }

}
