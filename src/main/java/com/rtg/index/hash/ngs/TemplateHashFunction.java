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

/**
 * Common functions to both <code>ReadHashFunction</code> and <code>TemplateHashFunction</code>.
 */
public interface TemplateHashFunction extends CommonHashFunction {

  /**
   * Compute fast score for the specified read against the current state of the hash function.
   * @param readId read identifier.
   * @return the score.
   */
  int fastScore(final int readId);

  /**
   * Compute indel score for the specified read against the current state of the hash function.
   * @param readId read identifier.
   * @return the score.
   */
  int indelScore(final int readId);

  /**
   * Set parameters for current template sequence..
   * @param name of the template.
   * @param length  of the template.
   */
  void templateSet(final long name, final int length);

  /**
   * Process all windows for the specified position in the template. This method processes forward and reverse frames.
   * @param endPosition last position in the current template window.
   * @throws IOException If an I/O error occurs
   */
  void templateBidirectional(final int endPosition) throws IOException;

  /**
   * Process all windows for the specified position in the template. This method processes forward frame only.
   * @param endPosition last position in the current template window.
   * @throws IOException If an I/O error occurs
   */
  void templateForward(final int endPosition) throws IOException;

  /**
   * Process all windows for the specified position in the template. This method processes reverse frame only.
   * @param endPosition last position in the current template window.
   * @throws IOException If an I/O error occurs
   */
  void templateReverse(final int endPosition) throws IOException;

  /**
   * Called when all positions in a sequence have been processed.
   * @throws IOException If an I/O error occurs
   */
  void endSequence() throws IOException;

    /**
   * Log statistics about performance.
   */
  void logStatistics();
}

