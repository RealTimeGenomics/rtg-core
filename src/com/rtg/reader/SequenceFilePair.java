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
package com.rtg.reader;

import java.io.IOException;

/**
 * internal type used by <code>SdfWriter</code> to manage writing sequence data
 */
interface SequenceFilePair {

  /**
   * Declares the beginning of a sequence.
   * @return false if operation was not performed due to not being enough space.
   * @throws IOException if an IO error occurs.
   */
  boolean markNextSequence() throws IOException;

  /**
   * Writes sequence bytes to disk.
   * @param data byte array to write
   * @param offset position in the array to start from (zero based)
   * @param length amount of data to write
   * @return true if succeeded, false if space limit hit.
   * @throws IOException if an I/O error occurs.
   */
  boolean write(byte[] data, int offset, int length) throws IOException;

  /**
   * Writes quality bytes to disk.
   * @param qual byte array to write
   * @param offset position in the array to start from (zero based)
   * @param length amount of data to write
   * @return true if succeeded, false if space limit hit.
   * @throws IOException if an I/O error occurs.
   */
  boolean writeQuality(byte[] qual, int offset, int length) throws IOException;

  /**
   * @return original length of data written
   */
  long valuesWritten();

  /**
   * @return the number of bytes available to write
   */
  long bytesFree();

  /**
   * Called to indicate have finished writing the final sequence.
   * @throws IOException if an I/O error occurs
   */
  void lastSequence() throws IOException;

  /**
   * Closes OutputStreams
   * @throws IOException if an I/O error occurs.
   */
  void close() throws IOException;

}
