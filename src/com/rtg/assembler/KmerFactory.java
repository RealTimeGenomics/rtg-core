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

import com.rtg.assembler.graph.Contig;

/**
 */
public interface KmerFactory {
  /**
   * Construct a Kmer
   *
   * @param kmer the bases this Kmer represents
   * @param start lowest base to include in the kmer
   * @param end highest base to include in the kmer
   * @return a new Kmer object for given string
   */
  Kmer make(byte[] kmer, int start, int end);
  /**
   * Construct a Kmer
   *
   * @param contig the bases this Kmer represents
   * @param start lowest base to include in the kmer
   * @param end highest base to include in the kmer
   * @return a new Kmer object for given string
   */
  Kmer make(Contig contig, int start, int end);
}
