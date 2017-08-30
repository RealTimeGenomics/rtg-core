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

package com.rtg.scheduler;

import com.rtg.util.MathUtils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class LookAhead extends IntegralAbstract {

  private final int mLookAhead;

  private final int mLength;

  private final int mMask;

  private final int[] mCounts;

  private final int mDelta;

  private int mCurrent = -1;

  private int mTotal = 0;


  @Override
  public synchronized boolean globalIntegrity() {
    integrity();
    if (mCurrent >= 0) {
      int total = 0;
      for (int i = mCurrent; i <= mCurrent + mLookAhead + mDelta; ++i) {
        final int cnt = mCounts[i & mMask];
        Exam.assertTrue(cnt >= 0);
        total += cnt;
      }
      Exam.assertEquals(total, mTotal);

      final int x = mCurrent & mMask;
      final int y = (mCurrent + mLookAhead + mDelta + 1) & mMask;
      if (x < y) {
        check(0, x);
        check(y, mLength);
      } else {
        check(y, x);
      }
    } else {
      check(0, mLength);
    }
    return true;
  }

  private void check(final int lo, final int hi) {
    for (int i = lo; i < hi; ++i) {
      Exam.assertEquals(0, mCounts[i & mMask]);
    }
  }


  @Override
  public synchronized boolean integrity() {
    Exam.assertTrue(mLookAhead >= 0);
    Exam.assertTrue(mLength > mLookAhead + 1);
    Exam.assertEquals(mMask, mLength - 1);
    Exam.assertTrue(mCurrent >= -1);
    Exam.assertTrue(mLookAhead + ":" + mTotal, mLookAhead + 1 >= mTotal && mTotal >= 0);
    if (mCurrent == -1) {
      Exam.assertEquals(0, mTotal);
    } else {
      Exam.assertTrue(mTotal >= 1);
    }
    if (mTotal == 0) {
      Exam.assertEquals(-1, mCurrent);
    }
    return true;
  }

  /**
   * @param lookAhead distance allowed to look ahead.
   * @param delta the distance from a task to the maximum of the tasks it is dependent on.
   */
  protected LookAhead(final int lookAhead, final int delta) {
    //System.err.println("LookAhead " + lookAhead);
    assert lookAhead >= 0;
    mLookAhead = lookAhead;
    mDelta = delta;
    mLength = 1 << MathUtils.ceilPowerOf2Bits(lookAhead + 1 + delta);
    //System.err.println(mLookAhead + ":" + mLength);
    mMask = mLength - 1;
    mCounts = new int[mLength];
  }

  /**
   * @return the total number of active times.
   */
  public synchronized int total() {
    return mTotal;
  }

  /**
   * @param t time to be checked.
   * @param delta additional distance forward that is ok.
   * @return true iff it is ok to schedule something at time t (no active entry closer than look ahead - no active entry <code>x &lt; t - lookAhead</code>).
   */
  public synchronized boolean ok(final int t, final int delta) {
    assert t >= 0;
    final boolean res;
    if (mCurrent == -1) {
      res = true;
    } else {
      res = mCurrent <= t && t <= mCurrent + mLookAhead + delta;
    }
    //System.err.println("LookAhead ok " + t + " " + res);
    return res;
  }

  /**
   * Mark time t as active.
   * @param t time to be marked.
   */
  public synchronized void increment(final int t) {
    //System.err.println("LookAhead increment " + t + " mCurrent=" + mCurrent);
    assert t >= 0;
    if (t < mCurrent || !ok(t, mDelta)) {
      throw new RuntimeException("Time invalid: current=" + mCurrent + " t=" + t + " lookAhead=" + mLookAhead + " mDelta=" + mDelta);
    }
    final int x = t & mMask;
    assert mCounts[x] >= 0;
    mCounts[x]++;
    ++mTotal;
    if (mCurrent == -1) {
      assert mTotal == 1;
      mCurrent = t;
    }
  }

  /**
   * Mark time t as no longer active.
   * @param t time to be marked.
   */
  synchronized void decrement(final int t) {
    //System.err.println("LookAhead decrement " + t + " mCurrent=" + mCurrent);
    assert t >= 0;
    final int x = t & mMask;
    assert mCounts[x] >= 1;
    mCounts[x]--;
    --mTotal;
    assert mCounts[x] >= 0;
    assert mTotal >= 0;
    if (mCurrent == t) {
      if (mTotal == 0) {
        mCurrent = -1;
      } else {
        while (true) {
          final int y = mCurrent & mMask;
          if (mCounts[y] > 0) {
            break;
          }
          ++mCurrent;
        }
      }
    }
  }

  /**
   * @return the look ahead.
   */
  public int lookAhead() {
    return mLookAhead;
  }

  /**
   * @return the delta.
   */
  public int delta() {
    return mDelta;
  }


}
