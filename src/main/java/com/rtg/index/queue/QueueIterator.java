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
