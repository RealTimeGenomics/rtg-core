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

package com.rtg.variant.bayes.multisample.family;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Provide unique identifiers in an unbroken sequence for a set of integer values.
 */
public final class UniqueId extends IntegralAbstract {

  private final int[] mIds;

  // -1 reserved externally for unseen but internally is 0
  private int mLastId = 0;

  /**
   * @param length the numbers which are to be identified range over <code>[0..length)</code>.
   */
  public UniqueId(final int length) {
    assert length >= 1;
    mIds = new int[length];
  }

  /**
   * Add another value to the list with unique identifiers.
   * @param i the value to be identified.
   * @return the (possibly new) identifier for <code>i</code> (&ge; 0).
   */
  public int addId(final int i) {
    final int id = mIds[i];
    if (id > 0) {
      return id - 1;
    }
    //increment before returning so that -1 reserved externally for unseen but internally is 0
    ++mLastId;
    mIds[i] = mLastId;
    return mLastId - 1;
  }

  /**
   * If i has been seen before return its identifier otherwise -1.
   * @param i the value to be identified.
   * @return the identifier for <code>i</code> (&ge; -1).
   */
  int id(final int i) {
    return mIds[i] - 1;
  }

  /**
   * @return the number of different unique values seen so far. The identifiers returned so far will all
   * be &lt; <code>numberIdsSoFar</code> and the next new value will return <code>numberIdsSoFar</code>.
   */
  int numberIdsSoFar() {
    return mLastId;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (final int id : mIds) {
      Exam.assertTrue(0 <= id && id <= mIds.length);
    }
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mIds.length >= 1);
    Exam.assertTrue(0 <= mLastId && mLastId <= mIds.length);
    return true;
  }

}
