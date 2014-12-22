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
 * Common interface between the array classes in the <code>longindex</code>,
 * <code>intindex</code> and <code>shortindex</code> packages.
 *
 */
public interface ExtensibleIndex extends CommonIndex {

  /**
   * Extend length by one and set newly created location to
   * value.
   * @param value the value as a long, regardless of the underlying type
   */
  void append(long value);

  /**
   * Extend length by one and set newly created location to
   * value.
   * Sign bits are preserved across <code>appendSigned</code> and <code>getSigned</code> calls.
   * @param value the value as a long, regardless of the underlying type
   */
  void appendSigned(long value);

  /**
   * Allocate an additional increment entries.
   * @param increment minimum number of new entries to be allocated.
   * @return initial position of the start of the newly allocated entries.
   */
  long extendBy(long increment);

  /**
   * Extend the array to support at least the specified size.
   * @param size new size
   */
  void extendTo(long size);

  /**
   * Reduce the length of the index to length and reduce the underlying memory used as far as possible.
   * The new length must be &le; the current length.
   * @param length new length of index.
   */
  void trim(long length);

  /**
   * Reduce the underlying memory used as far as possible to accomodate the current length.
   */
  void trim();

  /**
   * Gets the value at the specified position.
   * Sign bits are preserved across <code>setsigned</code> and <code>getSigned</code> calls.
   * This is not true of set and get calls where the low order bit pattern is preserved.
   * @param offset the position in this index
   * @return the value as a long, regardless of the underlying type
   */
  long getSigned(long offset);

  /**
   * Sets the value at the specified position.
   * Checks that the value is compatible with the underlying storage type.
   * In particular if the underlying storage type is n bits then a positive value
   * must have bit positions (n+1)..64 (inclusive, counting from 1) all 0.
   * If it is negative then it must have positions (n+1)..64 (inclusive, counting from 1) all 1.
   * This ensures that information is not lost.
   * Sign bits are preserved across <code>setsigned</code> and <code>getSigned</code> calls.
   * This is not true of set and get calls where the low order bit pattern is preserved.
   *
   * @param offset the position in this index
   * @param value the value as a long, regardless of the underlying type
   */
  void setSigned(long offset, long value);

}
