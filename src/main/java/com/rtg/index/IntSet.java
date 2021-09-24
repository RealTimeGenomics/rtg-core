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
package com.rtg.index;

import java.io.IOException;

import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public abstract class IntSet extends IntegralAbstract {

  private final IntSetCaller mCaller;

  /**
   * @param caller used to do calls for each unique value.
   */
  public IntSet(final IntSetCaller caller) {
    mCaller = caller;
  }

  /**
   * Add a new value to the set.
   * @param v the value to be added.
   * @throws IOException If an I/O error occurs
   */
  public abstract void add(final int v) throws IOException;

  /**
   * Call <code>call</code> for each unique value in the set.
   * Also remove all entries from the set.
   * @throws IOException If an I/O error occurs
   */
  public abstract void iterateClear() throws IOException;

  /**
   * Call <code>call</code> for each unique value in the set.
   * Also remove all entries from the set.
   * @throws IOException If an I/O error occurs
   */
  public abstract void iterateClearAll() throws IOException;

  /**
   * Called if the set overflows (exceeds capacity) or for each unique value
   * when <code>iterateClear</code> is called.
   * @param v value stored in set.
   * @throws IOException If an I/O error occurs
   */
  public void call(final int v) throws IOException {
    mCaller.call(v);
  }

}
