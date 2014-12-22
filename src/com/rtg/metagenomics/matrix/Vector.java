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
package com.rtg.metagenomics.matrix;


/**
 * One dimensional real valued vector.
 */
public interface Vector {

  /**
   * Get the <code>i'th</code> value.
   * @param i  index.
   * @return the value.
   */
  double get(int i);

  /**
   * Set the <code>i'th</code> value.
   * @param i index.
   * @param v value.
   */
  void set(int i, double v);

  /**
   * Increment the <code>i'th</code> value.
   * @param i index.
   * @param v value.
   */
  void incr(int i, double v);

  /**
   * Get the length of the vector.
   * @return the dimension (&ge; 0).
   */
  int dimension();
}
