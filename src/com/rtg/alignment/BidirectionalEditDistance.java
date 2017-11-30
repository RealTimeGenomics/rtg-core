/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.alignment;

/**
 * Interface for calculating an edit distance measure between two strings
 * in either forward or reverse directions
 */
public interface BidirectionalEditDistance {

  /**
   * Calculate a global alignment between an array and a specified portion
   * of another array. Higher weights, cause greater penalty for a particular
   * type of operation. There is an initial penalty for starting an indel
   * run, thereafter they suffer only ordinary penalty.
   * Ns added on either side of the template suffer an ordinary penalty.
   * @param read the read sequence
   * @param rlen length of read (assumed that given read is at least this long)
   * @param template the template sequence
   * @param zeroBasedStart the start position
   * @param rc true if a reverse complement alignment should be done
   * @param maxScore maximum score for an alignment to allow early termination
   * @param maxShift maximum allowed start position alignment shift
   * @param cgLeft true if this is the left arm (used only in CG case so far)
   * @return packed alignment information
   */
  int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, boolean rc, int maxScore, int maxShift, boolean cgLeft);

  /**
   * Log statistics from the alignment run.
   */
  void logStats();

}
