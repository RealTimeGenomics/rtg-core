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
