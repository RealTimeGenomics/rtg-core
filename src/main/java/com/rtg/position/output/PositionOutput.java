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
