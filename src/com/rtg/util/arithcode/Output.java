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

package com.rtg.util.arithcode;


/**
 */
public interface Output {

  /**
   * Get a key which is passed into Input to retrieve a particular block.
   * @return the key associated with the next block.
   */
  long endBlock();

  /**
   * Closes underlying output stream after filling to a byte
   * boundary with <code>0</code> bits.
   */
  void close();

  /**
   * Writes the single specified bit to the underlying output stream,
   * <code>1</code> for <code>true</code> and <code>0</code> for <code>false</code>.
   * @param bit Value to write.
   */
  void writeBit(boolean bit);

  /**
   * Writes a single <code>true</code> (<code>1</code>) bit.
   */
  void writeBitTrue();

  /**
   * Writes a single <code>false</code> (<code>0</code>) bit.
   */
  void writeBitFalse();

}
