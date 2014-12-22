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

package com.rtg.assembler;

import java.util.Iterator;

/**
 */
public interface DeBruijnGraph extends Iterable<Kmer> {
  /**
   * Add the provided Kmer to the graph or increase it's count
   * @param kmer the new kmer
   */
  //void add(Kmer kmer);

  //void add(KmerIterableFactory factory);

  /**
   * Call after adding all hashes to perform any compacting of finalisation of the graph before starting to use it
   */
  //void freeze();

  /**
   * Retrieve the number of times the kmer occurred in the input data
   * @param k the kmer to look up
   * @return  the number of times the kmer occurred in the input data
   */
  int frequency(Kmer k);

  /**
   * Should default to 1.
   * @param goodThreshold ignore kmer values which had less than this many occurrences in the input data
   */
  void setThreshold(int goodThreshold);

  /**
   * @param k the kmer value to update
   * @param built has the kmer been added to a contig
   */
  void setBuilt(Kmer k, boolean built);

  /**
   * @param k the kmer value to look up
   * @return true if the kmer is already part of a contig
   */
  boolean isBuilt(Kmer k);

  /**
   * @param k the kmer to look up
   * @return true if the kmer appears in the graph and is above threshold
   */
  boolean contains(Kmer k);

  /**
   * @return An iterator over all kmer in the graph that are above threshold
   */
  @Override
  Iterator<Kmer> iterator();

  /**
   * @return number of bytes used by the data structure. May neglect <code>O(1)</code> terms.
   */
  long bytes();
}
