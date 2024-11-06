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
package com.rtg.variant.match;

import com.rtg.util.Utils;
import com.rtg.variant.NoQualityException;

/**
 * Holds information about a match between a read and a reference
 * genome.
 */
public abstract class Match {

 /**
   * Get the base error information at position index.
   * This is the probability that the mapping of the corresponding nucleotide is incorrect
   * and is typically directly derived from the sequencer quality values.
   * @param index position to get quality from. Zero based.
   * @return the probability that the mapping is incorrect (may be 0.0 or 1.0).
   * @throws NoQualityException if no quality information is available
   *         (<code>hasQuality()</code> returns false).
   * @throws IndexOutOfBoundsException if index is negative or greater than or equal to <code>length()</code>.
   */
  public abstract double baseError(int index);

  /**
   * Get the average base error over the whole region.
   * This is the probability that the mapping of the corresponding nucleotides is incorrect
   * and is typically directly derived from the sequencer quality values.
   * @return the probability that the mapping is incorrect (may be 0.0 or 1.0).
   */
  public abstract double baseError();

  /**
   * Return a human readable version of the various quality log probabilities.
   * @return the string with quality information.
   */
  public String qualityString() {
    final StringBuilder sb = new StringBuilder();
    sb.append("[");
    for (int i = 0; i < length(); ++i) {
      sb.append(Utils.realFormat(-Math.log(baseError(i)) / Math.log(10.0), 2)).append(", ");
    }
    sb.append("]");
    return sb.toString();
  }

  /**
   * Get the read nucleotide at position index.
   *
   * @param index position to get nucleotide from. Zero based.
   * @return read nucleotide (range -1/N, 0/A, 1/C, 2/G, 3/T, 4/D).
   * @throws IndexOutOfBoundsException if index is negative or greater than or equal to <code>length()</code>.
   */
  public abstract int read(int index);

  /**
   * Get the probability that the overall mapping is incorrect.
   * @return the probability that the overall mapping is incorrect (can be 0.0 or 1.0).
   */
  public abstract double mapError();

  /**
   * Number of nucleotides (and quality values if they are present) in the read.
   * @return the number of nucleotide (zero or greater).
   */
  public abstract int length();

  /**
   * Return the read as a string representing the nucleotides as "" for a zero length string, and
   * A, C, G, T and D characters for individual nucleotides and deletions.
   * @return the read.
   */
  public abstract String readString();

  /**
   * Get the status of the match at its left end. If the left end is the end of the read then this is false.
   * @return left end status.
   */
  public abstract boolean isFixedLeft();

  /**
   * Get the status of the match at its right end. If the right end is the end of the read then this is false.
   * @return right end status.
   */
  public abstract boolean isFixedRight();


  /**
   * Calculate the overall quality score for this match. This averages the quality probabilities for each nucleotide and
   * also takes into account the overall quality of the map.
   * @return quality as a probability.
   */
  public double correction() {
    final double mq = mapError();
    return baseError() * (1.0 - mq) + mq;
  }

  /**
   * Include "~" as necessary to represent un-anchored strings.
   * @return the description
   */
  public String toString() {
    final String str = readString();
    return (isFixedLeft() ? "" : "~") + str + (isFixedRight() ? "" : "~");
  }

  /**
   * Return count of alleles that produce this allele from site specific priors.
   * @return count of alleles
   */
  public int alleleCount() {
    return 0;
  }
}
