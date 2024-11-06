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
package com.rtg.index;

import java.io.IOException;

/**
 * Extends the basic index to allow hashes longer than 64 bits (that is more than one long).
 *
 */
public interface IndexExtended extends Index {
  /**
   * Adds an extended hash code to the index and stores the associated value.
   * The id must be &gt;= 0
   * @param hash the hash key
   * @param value to be associated with the key
   * @throws IllegalStateException if index has been frozen.
   */
  void add(long[] hash, long value);

  /**
   * Search for the supplied hash code.
   * A call is done to <code>found(id)</code> for each hit where id is the associated
   * id.
   * @param hash the hash
   * @param finder the finder
   * @throws IllegalStateException if index has not been frozen.
   * @throws IOException if the finder produces such an exception.
   */
  void search(long[] hash, Finder finder) throws IOException;

  /**
   * Iterate over all entries in this index.
   * A call is done to <code>found(hash, value)</code> for each entry.
   * @param finder the finder
   * @throws IllegalStateException if index has not been frozen.
   */
  void scanAll(FinderHashValueExtended finder);

  /**
   * Determine whether the index contains the supplied hash code.
   * @param hash the hash
   * @return true iff the index contains the hash.
   * @throws IllegalStateException if index has not been frozen.
   */
  boolean contains(long[] hash);

  /**
   * Search for the supplied hash code and return the number of hits.
   * @param hash the hash
   * @return the number of hits for the hash (&gt;= 0).
   * @throws IllegalStateException if index has not been frozen.
   */
  int count(long[] hash);

  /**
   * Find the index of the first occurrence of the hash, or negative if the hash is not in the index.
   * @param hash the hash to find.
   * @return internal location of the hash (if &lt; 0 then not found).
   * @throws IllegalStateException if index has not been frozen.
   */
  long first(long[] hash);

  /**
   * Get the hash at the found location (see search).
   * @param found internal location of the hash (&gt; 0).
   * @return the hash code.
   */
  long[] getHashExtended(long found);
}


