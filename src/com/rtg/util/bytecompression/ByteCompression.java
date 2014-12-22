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
 * Interface for classes holding a compressed set of
 * byte arrays.
 */
public interface ByteCompression {

  /**
   * Method to add a byte array to the set.
   * This is <b>not</b> safe for access from multiple threads.
   * @param buffer contains the bytes to be compressed and added.
   * @param offset the offset into the buffer of the start.
   * @param length the number of bytes from the buffer to use.
   */
  void add(byte[] buffer, int offset, int length);

  /**
   * Method to get the byte array at the given index.
   * This is safe for access from multiple threads.
   * @param buffer a byte array of sufficient size to contain uncompressed data.
   * @param index the index of the byte array to fetch.
   * @param offset the offset into the uncompressed data to start from.
   * @param length the length of the uncompressed data to fetch.
   */
  void get(byte[] buffer, long index, int offset, int length);

  /**
   * Trim internal arrays to minimise total memory.
   * Prevents further adds.
   */
  void freeze();

  /**
   * Get the approximate number of bytes used.
   * @return the approximate number of bytes used.
   */
  long bytes();
}
