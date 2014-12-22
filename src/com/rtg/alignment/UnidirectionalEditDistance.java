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
package com.rtg.alignment;

/**
 * Interface for calculating an edit distance measure between two strings
 * in the forward direction only.
 */
public interface UnidirectionalEditDistance {

  /**
   * Calculate a global alignment between an array and a specified portion
   * of another array. Higher weights, cause greater penalty for a particular
   * type of operation. There is an initial penalty for starting an indel
   * run, thereafter they suffer only ordinary penalty.
   * Ns added on either side of the template suffer an ordinary penalty.
   * Reverse complement unsupported by this method, must be converted first.
   * @param read the read sequence
   * @param rlen length of read (assumed that given read is at least this long)
   * @param template the template sequence
   * @param zeroBasedStart the start position
   * @param maxScore maximum score for an alignment to allow early termination
   * @param maxShift maximum allowed alignment shift
   * @param cgLeft true if this is the left arm (used only in CG case so far)
   * @return packed alignment information
   */
  int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft);

  /**
   * Calculate an edit distance for a portion of a read and template where only
   * the start position is known.
   * @param read the read sequence
   * @param readStartPos the start position along the read (zero based, inclusive)
   * @param readEndPos the end position along the read (zero based, exclusive)
   * @param template the template sequence
   * @param templateStartPos the start position along the template (zero based, inclusive)
   * @param maxScore maximum score for an alignment to allow early termination
   * @param maxShift maximum allowed alignment shift
   * @return packed alignment information
   */
  int[] calculateEditDistanceFixedStart(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos, int maxScore, int maxShift);

  /**
   * Calculate an edit distance for a portion of a read and template where only
   * the end position is known. May shift the start position.
   * @param read the read sequence
   * @param readStartPos the start position along the read (zero based, inclusive)
   * @param readEndPos the end position along the read (zero based, exclusive)
   * @param template the template sequence
   * @param templateExpectedStartPos the expected start position along the template (zero based, inclusive)
   * @param templateEndPos the end position along the template (zero based, exclusive)
   * @param maxScore maximum score for an alignment to allow early termination
   * @param maxShift maximum allowed start position alignment shift
   * @return packed alignment information
   */
  int[] calculateEditDistanceFixedEnd(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateExpectedStartPos, int templateEndPos, int maxScore, int maxShift);

  /**
   * Calculate an edit distance for a bounded portion of a read and template where
   * both start and end position are known.
   * @param read the read sequence
   * @param readStartPos the start position along the read (zero based, inclusive)
   * @param readEndPos the end position along the read (zero based, exclusive)
   * @param template the template sequence
   * @param templateStartPos the start position along the template (zero based, inclusive)
   * @param templateEndPos the end position along the template (zero based, exclusive)
   * @param maxScore maximum score for an alignment to allow early termination
   * @param maxShift maximum allowed alignment shift
   * @return packed alignment information
   */
  int[] calculateEditDistanceFixedBoth(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos, int templateEndPos, int maxScore, int maxShift);

  /**
   * Log statistics from the alignment run.
   */
  void logStats();

}
