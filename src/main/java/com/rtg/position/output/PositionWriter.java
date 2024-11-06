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
