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
package com.rtg.util.array;


/**
 */
public interface CommonIndex {

  /**
   * Swap the values at the two specified locations.
   *
   * @param index1 the first index to be swapped
   * @param index2 the second index to be swapped
   */
  void swap(final long index1, final long index2);

  /**
   * Note: the length returned by this function should be precisely controlled
   * by the user, as such it should not reflect the amount of memory allocated
   * internally. Instead it should reflect the number of elements requested.
   * @return the length of the array
   */
  long length();

  /**
   * The length of this index in bytes.
   *
   * @return the number of bytes
   */
  long bytes();

  /**
   * Gets the value at the specified position.
   * The exact bit pattern as stored in <code>set</code> is returned.
   * In particular stored values which have their high order (sign) bit
   * set in the underlying implementations will not be returned
   * with this sign bit extended to the full long.
   *
   * @param offset the position in this index
   * @return the value as a long, regardless of the underlying type
   */
  long get(final long offset);

  /**
   * Sets the value at the specified position.
   * Checks that the value is compatible with the underlying storage type.
   * In particular if the underlying storage type is n bits then the value
   * must have bit positions (n+1)..64 (inclusive, counting from 1) all 0.
   * This ensures that information is not lost.
   *
   * @param offset the position in this index
   * @param value the value as a long, regardless of the underlying type
   */
  void set(final long offset, final long value);

  /**
   * Interrogate whether this index implementation is safe from word tearing for get and set operations.
   * Specifically is it safe for different threads to update adjacent elements without risking concurrent update errors.
   * @see <a href="http://docs.oracle.com/javase/specs/jls/se7/html/jls-17.html#jls-17.6">JLS on word tearing</a>
   * @return true if this implementation is not subject to word tearing.
   */
  boolean safeFromWordTearing();
}
