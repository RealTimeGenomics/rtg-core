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
 */
public interface Add {

  /**
   * Adds a hash code  to the index and stores the associated id.
   * The id must be &gt;= 0
   * @param hash the hash key
   * @param id to be associated with the key (might be a sequence id or a position)
   * @throws IllegalStateException if index has been frozen.
   */
  void add(long hash, long id) throws IllegalStateException;

}
