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

package com.rtg.variant.bayes.multisample.multithread;

import java.util.LinkedList;
import java.util.List;

import com.rtg.scheduler.EventList;
import com.rtg.scheduler.LookAhead;
import com.rtg.scheduler.enumtime.EnumTimeId;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Circular array for time with list for the enum values (should be 0, 1 or 2 entries).
 * @param <C> the type of the pair.
 */
public class EventListMultiSample<C extends EnumTimeId<?>> extends IntegralAbstract implements EventList<C> {

  private List<C>[] getList(int length) {
    @SuppressWarnings("unchecked")
    final List<C>[] list = (List<C>[]) new List<?>[length];
    for (int i = 0; i < list.length; ++i) {
      list[i] = new LinkedList<>();
    }
    return list;
  }


  private int mLength;

  private int mFirst = 0; //0 based inclusive

  private int mLast = 0; //0 based exclusive

  private List<C>[] mItems;

  /**
   * Comment.
   */
  public EventListMultiSample() {
    mLength = 1;
    mItems = getList(mLength);
  }

  @Override
  public C next(final LookAhead lookAhead) {
    int first = mFirst % mLength;
    while (mItems[first].isEmpty() && mFirst != mLast) {
      ++mFirst;
      ++first;
      if (first == mLength) {
        first = 0;
      }
    }
    if (mFirst == mLast || !lookAhead.ok(mFirst, 0)) {
      return null;
    }
    final C res = first(mItems[first]);
    assert res != null;
    return res;
  }

  /**
   * Get the first item in the set and remove it.
   * @param list to be modified.
   * @return the first item in the set.
   */
  C first(List<C> list) {
    return list.remove(0);
  }

  @Override
  public void add(C id) {
    final int time = id.time();
    final int last;
    final int first;
    if (mLast == mFirst) {
      last = time + 1;
      first = time;
    } else if (time >= mLast) {
      last = time + 1;
      first = mFirst;
    } else if (time < mFirst) {
      last = mLast;
      first = time;
    } else {
      last = mLast;
      first = mFirst;
    }
    //System.err.println(" first=" + first + " last=" + last + " mLength=" + mLength);
    final int length = last - first;
    if (length >= mLength) {
      //oops buffer not big enough
      resize(last - first);
    }
    mLast = last;
    mFirst = first;
    final int ti = time % mLength;
    mItems[ti].add(id);
  }

  private void resize(final int target) {
    final int newLength = newLength(mLength, target);
    //System.err.println("newLength=" + newLength);
    final List<C>[] newItems = getList(newLength);
    for (int i = mFirst; i < mLast; ++i) {
      final int x = i % mLength;
      final int y = i % newLength;
      newItems[y] = mItems[x];
    }
    mItems = newItems;
    mLength = newLength;
  }

  static int newLength(final int length, final int target) {
    int newLength = length;
    while (newLength <= target) {
      newLength = (newLength * 3 + 1) / 2;
    }
    return newLength;
  }

  @Override
  public boolean contains(C id) {
    final int time = id.time();
    if (time < mFirst || time >= mLast) {
      return false;
    }
    final int ti = time % mLength;
    return mItems[ti].contains(id);
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (int i = mFirst; i < mLast; ++i) {
      for (final C id : mItems[i % mLength]) {
        Exam.assertEquals(i, id.time());
      }
    }
    final int fi = mFirst % mLength;
    final int la = mLast % mLength;
    if (mFirst == mLast) {
      for (int i = 0; i < mLength; ++i) {
        check(i);
      }
    } else {
      if (fi < la) {
        for (int i = 0; i < fi; ++i) {
          check(i);
        }
        for (int i = la; i < mLength; ++i) {
          check(i);
        }
      } else {
        for (int i = la; i < fi; ++i) {
          check(i);
        }
      }
    }
    return true;
  }

  private void check(int i) {
    Exam.assertEquals(0, mItems[i].size());
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mFirst <= mLast);
    Exam.assertTrue(mLast - mFirst <= mLength);
    Exam.assertEquals(mLength, mItems.length);
    return true;
  }

}
