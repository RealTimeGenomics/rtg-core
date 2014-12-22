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

