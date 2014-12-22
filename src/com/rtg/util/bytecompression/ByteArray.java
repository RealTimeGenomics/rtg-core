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
package com.rtg.util.bytecompression;

/**
 */
public abstract class ByteArray {
  /**
   * Get a single byte.
   * If multiple bytes are required, the array version of get should be used
   * rather than this method, since it is faster.
   *
   * @param offset which byte to get.
   * @return a byte
   */
  public abstract byte get(long offset);

  /**
   * Reads <code>count</code> bytes, starting at <code>offset</code>.
   * @param dest the array to copy into, starting from position 0.
   * @param offset the position to start reading from.
   * @param count how many bytes to cover.
   */
  public abstract void get(final byte[] dest, final long offset, final int count);

  /**
   * Set a single byte.
   * If multiple bytes must be set, the array version of set should be used
   * rather than this method, since it is faster.
   *
   * @param offset which byte to set.
   * @param value the new value to put into the array.
   */
  public abstract void set(long offset, byte value);

  /**
   * Writes <code>buffer[0 .. count-1]</code> into the byte array, starting at <code>offset</code>.
   * @param offset the position to copy to.
   * @param buffer the bytes to copy.
   * @param count how many bytes to copy.
   */
  public abstract void set(final long offset, final byte[] buffer, final int count);

  /**
   * Writes <code>buffer[0 .. count-1]</code> into the byte array, starting at <code>offset</code>.
   * @param offset the position to copy to.
   * @param buffer the bytes to copy.
   * @param bOffset offset into buffer to start copying from
   * @param count how many bytes to copy.
   */
  public abstract void set(final long offset, final byte[] buffer, final int bOffset, final int count);

  /** @return the total number of values that can be stored in this array. */
  public abstract long length();

  /** @return the total amount of memory used by this array. */
  public long bytes() {
    return length();  // by default
  }

  /**
   * A factory method that allocates a suitable subclass of ByteArray.
   *
   * @param size the required length
   * @return a subclass of ByteArray
   */
  public static ByteArray allocate(final long size) {
    if (size <= Integer.MAX_VALUE) {
      return new SingleByteArray((int) size);
    } else {
      return new MultiByteArray(size);
    }
  }
}
