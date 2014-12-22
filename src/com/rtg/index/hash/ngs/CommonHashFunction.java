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


/**
 * Common functions to both <code>ReadHashFunction</code> and <code>TemplateHashFunction</code>.
 */
public interface CommonHashFunction {

  /**
   * Append one more code to the end of the window.
   * @param code to be appended.
   */
  void hashStep(final byte code);

  /**
   * Hash step to allow internal registers to shift without adding a new code.
   */
  void hashStep();

  /**
   * Reset as if no <code>hashStep</code> calls had been done.
   */
  void reset();

  /**
   * The total number of windows.
   * @return the number of windows.
   */
  int numberWindows();

  /**
   * Get the total window size in nucleotides.
   * In some split windows this may be split up into smaller chunks.
   * @return total window size.
   */
  int windowSize();

  /**
   * Returns the length of reads being processed.
   * @return the read length.
   */
  int readLength();
}

