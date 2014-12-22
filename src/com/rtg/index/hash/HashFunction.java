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
package com.rtg.index.hash;


/**
 * Abstract class for hashing sequences of codes.
 * Be careful some of the implementations need to call reset()
 * at the beginning of each calculation others can be used incrementally
 * (only call reset() when unknown codes are seen).
 *
 */
public interface HashFunction {

  /**
   * Reset the Hash function to an initial state.
   */
  void reset();

  /**
   * Make one step in the hash calculation.
   * After window size calls to this routine the returned code will be a valid
   * hash code.
   *
   * @param code the next code (intended to be derived from the ordinal or individual residues).
   * @return the current hash value.
   */
  long hashStep(final byte code);

  /**
   * @return number of codes in a window
   */
  int getWindowSize();
}

