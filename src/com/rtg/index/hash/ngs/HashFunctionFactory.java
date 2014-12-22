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
 * Creates new copies of <code>SplitHashLoop</code> when requested.
 * Also supplies the number of windows used (needed so that the indexes can be created before
 * the actual <code>SplitHashLoop</code> is available.
 */
public interface HashFunctionFactory {

  /**
   * Get the number of windows that will be used by the <code>SplitHashLoop</code>.
   * @return the number of windows that will be used by the <code>SplitHashLoop</code>.
   */
  int numberWindows();

  /**
   * Get the number of bits recorded in each hash window.
   * @return the number of bits recorded in each hash window.
   */
  int hashBits();

  /**
   * Get the number of bits needed to uniquely represent a window (the hash may lose information).
   * @return the number of bits needed to uniquely represent a window.
   */
  int windowBits();

  /**
   * Get the window size
   * @return the size
   */
  int windowSize();

  /**
   * Create a new instance of a <code>HashFunction</code>.
   * @param readCall will do processing of read sequences (during build).
   * @param templateCall will do processing of template sequences (during search).
   * @return a new instance of a <code>HashFunction</code>.
   */
  NgsHashFunction create(ReadCall readCall, TemplateCall templateCall);
}

