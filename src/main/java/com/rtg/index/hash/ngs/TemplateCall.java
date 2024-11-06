/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

