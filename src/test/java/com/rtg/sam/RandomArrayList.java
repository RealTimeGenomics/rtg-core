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
