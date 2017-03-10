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
package com.rtg.position.output;

import java.io.IOException;

import com.rtg.reader.NamesInterface;

/**
 * Interface defining a surrogate object capable of printing a result.
 *
 */
public interface SurrogateRegion {

  /**
   * Write a result.
   *
   * @param out place where result should be written
   * @param subjectNames Table for subject names
   * @param queryNames Table for query names
   * @return true if output is actually written (useful for determining if a
   * result exceeding any score threshold)
   * @throws IOException if an IO error occurs
   */
  boolean write(final Appendable out, final NamesInterface subjectNames, final NamesInterface queryNames) throws IOException;

  /**
   * Write the header in the result
   * @param out place where header should be written
   * @throws IOException if an IO error occurs
   */
  void writeHeader(final Appendable out) throws IOException;

  /**
   * Returns <code>false</code> if this score doesn't pass the threshold test.
   * @return true for keep
   */
  boolean scoreAllowed();

  /**
   * Return the score for the underlying region.
   *
   * @return score
   */
  double score();

  /**
   * Return the query id.
   *
   * @return query id
   */
  int queryId();
}
