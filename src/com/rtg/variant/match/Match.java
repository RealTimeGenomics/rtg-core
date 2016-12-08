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
  public abstract double baseError(int index) throws NoQualityException, IndexOutOfBoundsException;

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
