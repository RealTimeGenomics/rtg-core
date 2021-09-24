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
package com.rtg.index.hash.ngs;

import java.io.IOException;

import com.rtg.launcher.HashingRegion;

/**
 * Used for he operations done by the <code>windowAll</code> method.
 */
public interface TemplateCall extends Cloneable {
  /**
   * Do processing for any hits from the specified hash code against the selected index.
   * @param endPosition position of the last residue of the current window.
   * @param hash generated for the window.
   * @param index of the window used.
   * @throws IOException If an I/O error occurs
   */
  void templateCall(int endPosition, long hash, int index) throws IOException;

  /**
   * Set the parameters for the template sequence currently being processed.
   * @param name of the sequence.
   * @param length of the sequence.
   */
  void set(long name, int length);

  /**
   * Specify the hash function available for use (for example for scoring).
   * @param hashFunction to be used.
   */
  void setHashFunction(NgsHashFunction hashFunction);

  /**
   * Called when all hits for a particular position have been processed.
   * @throws IOException If an I/O error occurs
   */
  void done() throws IOException;

  /**
   * Called when all hits for a particular sequence have been processed.
   * @throws IOException If an I/O error occurs
   */
  void endSequence() throws IOException;

  /**
   * Sets whether we are operating in the reverse frame.
   * @param reverse true iff we are in reverse frame.
   */
  void setReverse(boolean reverse);

  /**
   * Returns true if we are processing reverse frame.
   * @return true iff reverse frame.
   */
  boolean isReverse();

  /**
   * Create a clone of this TemplateCall that is suitable for running in an independent
   * thread. This method will ensure that any shared configuration is copied or that
   * thread-safe access is provided.
   *
   * @param region the sequence region beyond which results should be filtered.
   * @return a clone.
   * @throws IOException If an I/O error occurs
   */
  TemplateCall threadClone(HashingRegion region) throws IOException;

  /**
   * Close per thread resources
   * @throws IOException if an IO error occurs
   */
  void threadFinish() throws IOException;

  /**
   * Create a clone.
   * @return a clone.
   * @throws CloneNotSupportedException when error from <code>Object</code> level clone.
   */
  TemplateCall clone() throws CloneNotSupportedException;

  /**
   * Output to log accumulated statistics about performance.
   */
  void logStatistics();

}

