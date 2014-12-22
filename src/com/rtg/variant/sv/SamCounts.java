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

package com.rtg.variant.sv;

/**
 * Records counts of SAM records at particular coordinates.
 *
 */
public interface SamCounts {

  /**
   * Add a new record to the underlying set at position given by index.
   * @param index position used.
   */
  void increment(final int index);

  /**
   * Get the number of records at the index.
   * @param base position for window (in forward case this is at right end of
   * window).
   * @param index location being fetched (in forward case counting from left hand end of window).
   * @return count - 0 if no set at specified location
   */
  double count(int base, int index);

  /**
   * Get the number of records at all positions between <code>index</code> and <code>index2</code> (exclusive).
   * @param base position for window.
   * @param index first location being fetched (in forward case counting from left hand end of window).
   * @param index2 last location being fetched (exclusive).
   * @return count - 0 if no set at specified location
   */
  double count(int base, int index, int index2);

  /**
   * Get the sum of count * log(count) between <code>index</code> and <code>index2</code> (exclusive).
   * @param base position for window.
   * @param index first location being fetched (in forward case counting from left hand end of window).
   * @param index2 last location being fetched (exclusive).
   * @return count - 0 if no set at specified location
   */
  double sumLn(int base, int index, int index2);

  /**
   * Reset all sets at all indexes and ensure that the array
   * is able to accommodate length entries.
   * @param length to be accommodated from now on.
   * @param bufferSize size of active view
   */
  void reset(final int length, int bufferSize);

  /**
   * Sets the current minimum position that can be used
   * @param offset the current minimum position.
   */
  void flushTo(int offset);

  /**
   * The current effective length of the underlying index (used for checking index values etc.).
   * May be less than the actual backing array.
   * @return the effective length.
   */
  int length();

}
