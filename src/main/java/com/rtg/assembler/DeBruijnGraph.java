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

package com.rtg.assembler;

import java.util.Iterator;

/**
 * Definition of graph used in assembler.
 */
public interface DeBruijnGraph extends Iterable<Kmer> {

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
