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

package com.rtg.sam;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * @param <A> type of items in the array list.
 */
public class RandomArrayList<A> extends IntegralAbstract {

  private final ArrayList<A> mList = new ArrayList<>();

  private final Random mRandom;

  /**
   * @param seed for the psueudo-random number generator.
   */
  public RandomArrayList(final long seed) {
    mRandom = new Random(seed);
  }

  /**
   * @param seed for the pseudo-random number generator.
   * @param items to be added to the array list.
   */
  public RandomArrayList(final long seed, final A[] items) {
    this(seed);
    Collections.addAll(mList, items);
  }

  /**
   * Add a new item to the array list.
   * @param item to be added.
   */
  public void add(final A item) {
    mList.add(item);
  }

  /**
   * Remove and return a random entry from the array list.
   * @return the next random entry from the array list.
   */
  public A next() {
    if (mList.isEmpty()) {
      return null;
    }
    final int size = mList.size();
    final int nx = mRandom.nextInt(size);
    final A res = mList.get(nx);
    final A last = mList.remove(size - 1);
    if (nx != size - 1) {
      mList.set(nx, last);
    }
    return res;
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mList);
    Exam.assertNotNull(mRandom);
    return true;
  }

}
