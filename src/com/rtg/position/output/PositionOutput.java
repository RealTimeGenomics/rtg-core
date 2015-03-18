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
 */
public interface PositionOutput {

  /**
   * Called when there is a hit in the index.
   * @param seqId identifier of the sequence (from build).
   * @param posn position within the sequence (starts at 0, start of hit).
   * @throws IOException when writing output.
   */
  void hit(int seqId, int posn) throws IOException;

  /**
   * Set the current value of the position in the search sequence.
   * @param position the current value of the position in the search sequence.
   * Starts at 0.
   * @throws IOException when writing output.
   */
  void setPosition(int position) throws IOException;

  /**
   * Called when all the hits for a particular position have been finished.
   * @throws IOException if any underlying writing causes an error.
   */
  void endPosition() throws IOException;

  /**
   * Called when step from one sequence to next in search.
   * @param seqId query sequence id
   * @param length of the current query sequence.
   * @throws IOException if an IO error occurs
   */
  void nextSequence(int seqId, int length) throws IOException;

  /**
   * Called when step from one sequence to next in search.
   * @param seqId query sequence id
   * @param length of the current query sequence (that is, the entire sequence not just the part that is represented in <code>sequence</code>).
   * @param usedLength initial length of <code>sequence</code> that is actually used.
   * @param sequence the query sequence of nucleotides.
   * @throws IOException if an IO error occurs
   */
  void nextSequence(int seqId, int length, int usedLength, byte[] sequence) throws IOException;

  /**
   * Called when step from one frame to next in search.
   * @param frame framing of the search sequence.
   * @param seqId identifier of the search sequence.
   */
  void nextQuery(Frame frame, int seqId);

  /**
   * Called when all the positions for a particular query have been finished.
   * @throws IOException if any underlying writing causes an error.
   */
  void endQuery() throws IOException;

  /**
   * Called when all positions for a particular query/frame pair have been finished.
   * @throws IOException if any underlying writing causes an error.
   */
  void endQuerySequence() throws IOException;

  /**
   * Called at the end when all processing has been done.
   * @throws IOException if any underlying writing causes an error.
   */
  void endAll() throws IOException;

  /**
   * Get the total score accumulated during output.
   * Details of how it is computed vary between the different implementations.
   * @return total score.
   */
  double score();

  /**
   * Creates a clone for handling reverse frame calls in parallel. General contract
   * implies caller will not make repeat <code>nextQuery/endAll</code> etc calls if this is same
   * object, otherwise will.
   * @return the clone
   */
  PositionOutput reverseClone();
}
