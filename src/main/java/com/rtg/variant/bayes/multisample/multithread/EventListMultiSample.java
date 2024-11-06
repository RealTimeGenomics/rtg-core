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
    while (mFirst != mLast && mItems[first].isEmpty()) {
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
