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

import com.rtg.mode.Frame;

/**
 * Interface for classes taking output from <code>Position</code> and selecting
 * all or various subsets for output.
 *
 */
public interface PositionWriter {

  /**
   * End processing of query current sequence.
   *
   * @param queryId the query id of the sequence just finished
   * @throws IOException If an IO error occurs
   */
  void endQuery(final int queryId) throws IOException;

  /**
   * Called once hits for all queries have been found.
   */
  void endAll();

  /**
   * Write a result.
   *
   * @param region region to write
   * @param sequenceId sequence id
   * @param queryFrame frame of the query sequence
   * @param queryLength length of the query sequence
   * @param queryEffectiveLength number of words that can be fitted into the query sequence (adjusted for word length).
   * @throws IOException If an IO Error occurs
   */
  void write(final AbstractGappedRegion<?> region, final int sequenceId, final Frame queryFrame, final int queryLength, int queryEffectiveLength) throws IOException;

  /**
   * Return the total score of all written results.
   *
   * @return total score
   */
  double score();

  /**
   * Creates a clone for performing operations on a another frame in parallel
   * @return the clone
   */
  PositionWriter reverseClone();
}
