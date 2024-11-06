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
