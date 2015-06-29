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

package com.rtg.util.memo;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;


/**
 * Remember values of integer functions defined over all integer values (including
 * negative ones).
 * Assumes that the underlying function is well defined (doesn't throw any exceptions)
 * on a contiguous region of values.
 */
public class DoubleMemo extends IntegralAbstract implements DoubleFunction {

  // Care is taken to make this thread safe which includes storing state in the class Pair and updating it
  // in a single atomic operation. Also there is a short circuit that returns results in the common case
  // when the function value is available, this requires the state in an atom so that it can be coherently accessed.
  // Correctness of this assumes pointer updates are atomic and that the JIT cannot rearrange the order of operations
  // when executing a constructor and then assigning the result (that is all stores in the constructor occur before
  // the resultant object is itself stored).
  private static class Pair {
    private final int mLo;
    private final double[] mMemo;

    Pair(int lo, double[] memo) {
      mLo = lo;
      mMemo = memo;
    }
  }

  private static double[] makeArray(final int length) {
    return new double[length];
  }

  private final DoubleFunction mFunction;

  private Pair mSynchPair = null;

  /**
   * @param function underlying function being reflected in memo.
   */
  public DoubleMemo(final DoubleFunction function) {
    mFunction = function;
  }

  @Override
  public double fn(final int i) {
    final Pair memo = mSynchPair;
    if (memo != null) {
      //special case that delicately avoids synchronization.
      //memo may be an old version but that doesn't matter
      final int index = i - memo.mLo;
      if (index >= 0 && index < memo.mMemo.length) {
        return memo.mMemo[index];
      }
    }
    synchronized (this) {
      final int ix;
      if (mSynchPair == null) {
        final double[] memArray = makeArray(1);
        memArray[0] = mFunction.fn(i);
        mSynchPair = new Pair(i, memArray);
        ix = 0;
      } else {
        final int index = i - mSynchPair.mLo;
        if (index < 0) {
          final int newLength = mSynchPair.mMemo.length - index;
          final double[] newArray = makeArray(newLength);
          System.arraycopy(mSynchPair.mMemo, 0, newArray, -index, mSynchPair.mMemo.length);
          for (int j = 0; j < -index; j++) {
            newArray[j] = mFunction.fn(j + i);
          }
          mSynchPair = new Pair(i, newArray);
          ix = 0;
        } else if (index >= mSynchPair.mMemo.length) {
          final int newLength = index + 1;
          final double[] newArray = makeArray(newLength);
          System.arraycopy(mSynchPair.mMemo, 0, newArray, 0, mSynchPair.mMemo.length);
          for (int j = mSynchPair.mMemo.length; j <= index; j++) {
            newArray[j] = mFunction.fn(j + mSynchPair.mLo);
          }
          mSynchPair = new Pair(mSynchPair.mLo, newArray);
          ix = index;
        } else {
          //if we come through here then we have dodged a race condition
          ix = index;
        }
      }
      return mSynchPair.mMemo[ix];
    }
  }

  @Override
  public synchronized void toString(StringBuilder sb) {
    sb.append("IntMemo ");
    if (mSynchPair == null) {
      sb.append("empty");
    } else {
      sb.append("[").append(mSynchPair.mLo).append("..").append(mSynchPair.mLo + mSynchPair.mMemo.length - 1).append("]").append(LS);
      boolean first = true;
      for (final double fn : mSynchPair.mMemo) {
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
      for (int i = 0; i < mSynchPair.mMemo.length; i++) {
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
