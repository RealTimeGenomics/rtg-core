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

package com.rtg.util.memo;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;


/**
 * Remember values of integer functions defined over all integer values (including
 * negative ones).
 * Assumes that the underlying function is well defined (doesn't throw any exceptions)
 * on a contiguous region of values.
 * @param <T> type returned by the function.
 */
public class Memo<T> extends IntegralAbstract implements Function<T> {
  // Care is taken to make this thread safe which includes storing state in the class Pair and updating it
  // in a single atomic operation. Also there is a short circuit that returns results in the common case
  // when the function value is available, this requires the state in an atom so that it can be coherently accessed.
  // Correctness of this assumes pointer updates are atomic and that the JIT cannot rearrange the order of operations
  // when executing a constructor and then assigning the result (that is all stores in the constructor occur before
  // the resultant object is itself stored).
  private static class Pair {
    private final int mLo;
    private final Object[] mMemo;

    Pair(int lo, Object[] memo) {
      mLo = lo;
      mMemo = memo;
    }
  }

  private final Function<T> mFunction;

  private Pair mSynchPair = null;

  /**
   * @param function underlying function being reflected in memo.
   */
  public Memo(final Function<T> function) {
    mFunction = function;
  }

  @Override
  public T fn(final int i) {
    final Pair memo = mSynchPair;
    if (memo != null) {
      //special case that delicately avoids synchronization.
      //memo may be an old version but that doesn't matter
      final int index = i - memo.mLo;
      if (index >= 0 && index < memo.mMemo.length) {
        @SuppressWarnings("unchecked")
        final T ret = (T) memo.mMemo[index];
        return ret;
      }
    }
    synchronized (this) {
      final int ix;
      if (mSynchPair == null) {
        final Object[] memArray = new Object[1];
        memArray[0] = mFunction.fn(i);
        mSynchPair = new Pair(i, memArray);
        ix = 0;
      } else {
        final int index = i - mSynchPair.mLo;
        if (index < 0) {
          final int newLength = mSynchPair.mMemo.length - index;
          final Object[] newArray = new Object[newLength];
          System.arraycopy(mSynchPair.mMemo, 0, newArray, -index, mSynchPair.mMemo.length);
          for (int j = 0; j < -index; ++j) {
            newArray[j] = mFunction.fn(j + i);
          }
          mSynchPair = new Pair(i, newArray);
          ix = 0;
        } else if (index >= mSynchPair.mMemo.length) {
          final int newLength = index + 1;
          final Object[] newArray = new Object[newLength];
          System.arraycopy(mSynchPair.mMemo, 0, newArray, 0, mSynchPair.mMemo.length);
          for (int j = mSynchPair.mMemo.length; j <= index; ++j) {
            newArray[j] = mFunction.fn(j + mSynchPair.mLo);
          }
          mSynchPair = new Pair(mSynchPair.mLo, newArray);
          ix = index;
        } else {
          //if we come through here then we have dodged a race condition
          ix = index;
        }
      }
      @SuppressWarnings("unchecked")
      final T ret = (T) mSynchPair.mMemo[ix];
      return ret;
    }
  }

  @Override
  public synchronized void toString(StringBuilder sb) {
    sb.append("Memo ");
    if (mSynchPair == null) {
      sb.append("empty");
    } else {
      sb.append("[").append(mSynchPair.mLo).append("..").append(mSynchPair.mLo + mSynchPair.mMemo.length - 1).append("]").append(LS);
      boolean first = true;
      for (final Object fn : mSynchPair.mMemo) {
        if (!first) {
          sb.append(", ");
        }
        sb.append("").append(fn);
        first = false;
      }
    }
    sb.append(LS);
  }

  @Override
  public synchronized boolean globalIntegrity() {
    integrity();
    if (mSynchPair != null) {
      for (int i = 0; i < mSynchPair.mMemo.length; ++i) {
        final int j = i + mSynchPair.mLo;
        Exam.assertEquals("i", mFunction.fn(j), mSynchPair.mMemo[i]);
      }
    }
    return true;
  }

  @Override
  public synchronized boolean integrity() {
    Exam.assertNotNull(mFunction);
    return true;
  }

}
