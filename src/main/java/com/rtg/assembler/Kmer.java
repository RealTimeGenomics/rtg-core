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

/**
 */
interface Kmer extends Comparable<Kmer> {
  /**
   * Construct a new <code>Kmer</code> of the same type but representing one base right
   * @param nt the new base to append in the range A..T (1..4)
   * @return the kmer with <code>nt</code> appended and the first base dropped
   */
  Kmer successor(byte nt);

  /**
   * Construct a new <code>Kmer</code> of the same type but representing one base left
   * @param nt the new base to prepend in the range A..T (1..4)
   * @return the kmer with <code>nt</code> prepended and the last base dropped
   */
  Kmer predecessor(byte nt);

  /**
   * Consistently choose one of this Kmer or it's reverse complement.
   * If K is a kmer and R is it's reverse complement then:
   * <code>K.minimalKmer() == R.minimalKmer()</code> and
   * <code>K.minimalKmer() == K || K.minimalKmer() == R</code>
   * @return The 'standard' of this kmer and it's reverse complement
   */
  Kmer minimalKmer();

  /**
   *
   * @return a kmer representing the reverse complement of this one
   */
  Kmer reverse();

  /**
   * @return the length of this kmer in bases
   */
  int length();

  /**
   * @param pos a position within this kmer
   * @return the base at position <code>pos</code> encoded according to <code>DNA</code>
   */
  byte nt(int pos);
}
