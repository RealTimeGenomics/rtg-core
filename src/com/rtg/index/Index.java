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
  void freeze() throws IllegalStateException;

  /**
   * Search for the supplied hash code.
   * A call is done to <code>found(id)</code> for each hit where id is the associated
   * id.
   * @param hash the hash
   * @param finder the finder
   * @throws IllegalStateException if index has not been frozen.
   * @throws IOException if the finder produces such an exception.
   */
  void search(long hash, Finder finder) throws IOException, IllegalStateException;

  /**
   * Iterate over all entries in this index.
   * A call is done to <code>found(hash, value)</code> for each entry.
   * @param finder the finder
   * @throws IllegalStateException if index has not been frozen.
   * @throws IOException if the finder produces such an exception.
   */
  void scan(FinderHashValue finder) throws IOException, IllegalStateException;

  /**
   * Determine whether the index contains the supplied hash code.
   * @param hash the hash
   * @return true iff the hash is contained in the index
   * @throws IllegalStateException if index has not been frozen.
   */
  boolean contains(long hash) throws IllegalStateException;

  /**
   * Search for the supplied hash code and return the number of hits.
   * @param hash the hash
   * @return the number of hits for the hash (&gt;= 0).
   * @throws IllegalStateException if index has not been frozen.
   */
  int count(long hash) throws IllegalStateException;

  /**
   * Find the index of the first occurrence of the hash, or negative if the hash is not in the index.
   * @param hash the hash to find.
   * @return internal location of the hash (if &lt; 0 then not found).
   * @throws IllegalStateException if index has not been frozen.
   */
  long first(long hash) throws IllegalStateException;

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
  long numberEntries() throws IllegalStateException;

  /**
   * Get the number of hashes actually stored.
   * This is after duplicates and hashes which contain more than repeat
   * frequency threshold entries have been removed.
   * @return the number of hashes stored
   * @throws IllegalStateException if index is not frozen.
   */
  long numberHashes() throws IllegalStateException;

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
}


