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
package com.rtg.index;

/**
 * specifies the interface used for those filtering hashes during index compact stage
 */
public interface IndexFilterMethod {

  /**
   * Called once before using
   * @param index the index we will be filtering hashes from
   */
  void initialize(Index index);

  /**
   * @param hash the hash
   * @param numHits the number of occurences of this hash in the index
   * @return true if hash should be kept, false if it should be discarded
   */
  boolean keepHash(long hash, long numHits);

  /**
   * @return a copy of this filter
   */
  IndexFilterMethod threadClone();
}
