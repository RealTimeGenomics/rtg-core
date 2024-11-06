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
