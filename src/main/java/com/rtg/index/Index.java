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
import java.io.PrintStream;

/**
 * Interface for hash-based sequence indexes.
 *
 * The order that calls can be made is:
 * <code> (add* freezeIndex search* melt)* [add* freezeIndex search*]</code>
 *
 */
public interface Index extends Add {

  /**
   * Performs any operations that might be required to get the index
   * into a state where it can be queried. This method must be called
   * before any searching is to be done.
   * It is not possible to call add after this.
   * @throws IllegalStateException if index has already been frozen.
   */
  void freeze();

  /**
   * Search for the supplied hash code.
   * A call is done to <code>found(id)</code> for each hit where id is the associated
   * id.
   * @param hash the hash
   * @param finder the finder
   * @throws IllegalStateException if index has not been frozen.
   * @throws IOException if the finder produces such an exception.
   */
  void search(long hash, Finder finder) throws IOException;

  /**
   * Iterate over all entries in this index.
   * A call is done to <code>found(hash, value)</code> for each entry.
   * @param finder the finder
   * @throws IllegalStateException if index has not been frozen.
   */
  void scan(FinderHashValue finder);

  /**
   * Determine whether the index contains the supplied hash code.
   * @param hash the hash
   * @return true iff the hash is contained in the index
   * @throws IllegalStateException if index has not been frozen.
   */
  boolean contains(long hash);

  /**
   * Search for the supplied hash code and return the number of hits.
   * @param hash the hash
   * @return the number of hits for the hash (&gt;= 0).
   * @throws IllegalStateException if index has not been frozen.
   */
  int count(long hash);

  /**
   * Returns a human readable string reporting on
   * performance statistics.
   *
   * @return a human readable string reporting on memory usage and
   * performance statistics.
   */
  String perfString();

  /**
   * Returns a human readable string reporting on
   * memory usage statistics.
   *
   * @return a human readable string reporting on memory usage statistics.
   */
  String infoString();

  /**
   * Calculate the total memory used by the various parts of the index.
   * Only takes into account the parts that are determined by size -
   * constant terms are ignored.
   * @return the total memory used in bytes.
   */
  long bytes();

  /**
   * Get the number of distinct entries (hash, value pairs) stored.
   * This is after duplicates and hashes which contain more than repeat
   * frequency threshold entries have been removed.
   * @return the number of entries stored
   * @throws IllegalStateException if index is not frozen.
   */
  long numberEntries();

  /**
   * Get the number of hashes actually stored.
   * This is after duplicates and hashes which contain more than repeat
   * frequency threshold entries have been removed.
   * @return the number of hashes stored
   * @throws IllegalStateException if index is not frozen.
   */
  long numberHashes();

  /**
   * @return the number of hashes added before the freeze stage
   */
  long getInitialHashes();

  /**
   * Will write out hashes and values in a format identical for both overflow and no overflow cases.
   * Needed to investigate a bug.
   * @param out where to put output.
   */
  void dumpValues(final PrintStream out);

  /**
   * @return Number of occurences of the most frequent hash in the index
   */
  int maxHashCount();

  /**
   * Find the index of the first occurrence of the hash, or negative if the hash is not in the index.
   * @param hash the hash to find.
   * @return internal location of the hash (if &lt; 0 then not found).
   * @throws IllegalStateException if index has not been frozen.
   */
  long first(long hash);

  /**
   * Get the hash at the found location (see search).
   * @param found internal location of the hash (&gt; 0).
   * @return the hash code.
   */
  long getHash(long found);

  /**
   * Get the value associated with the hash at found location (see search). Note that search will only
   * return a single location even though multiple may be associated with the hash.
   * @param found internal location of value (&gt;= 0).
   * @return the value
   */
  long getValue(long found);


  /**
   * @return histogram of hash frequencies
   */
  SparseFrequencyHistogram getSparseFrequencyHistogram();
}


