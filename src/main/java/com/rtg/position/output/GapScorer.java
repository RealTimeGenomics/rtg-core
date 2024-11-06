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

import com.rtg.mode.Frame;


/**
 */
public interface GapScorer {

  /**
   * Provides query id.
   * @param queryId id of current query as from {@link com.rtg.reader.SequencesReader}
   * @param frame frame of query.
   */
  void setQueryId(int queryId, Frame frame);

  /**
   * Provides internal build id
   * @param internalId id of current build sequence encoded for frame
   */
  void setInternalBuildId(int internalId);

  /**
   * Return the score for the gap. Bigger = better
   * Note: in the reverse complement case, the read positions are reported as if
   * the read is reversed, and the template positions as if the template is forward.
   * @param buildStart the start of the region we wish to score on the build (incl)
   * @param buildEnd the end of the region we wish to score on the build (excl)
   * @param queryStart the start of the region we wish to score on the query (incl)
   * @param queryEnd the end of the region we wish to score on the query (excl)
   * @return the score
   */
  double score(final int buildStart, final int buildEnd, final int queryStart, final int queryEnd);

  /**
   * Return the score for the gaps at the ends. Bigger = better
   * Note: in the reverse complement case, the read positions are reported as if
   * the read is reversed, and the template positions as if the template is forward.
   * @param buildStart the start of the region we wish to score on the build (incl)
   * @param buildEnd the end of the region we wish to score on the build (excl)
   * @param queryStart the start of the region we wish to score on the query (incl)
   * @param queryEnd the end of the region we wish to score on the query (excl)
   * @return a maximum score
   */
  double scoreMax(final int buildStart, final int buildEnd, int queryStart, int queryEnd);

  /**
   * Find the maximum delta value that will give a non-zero score.
   * @return the maximum delta (likely to be &gt; 0 but this isn't guaranteed).
   */
  int maxDelta();

  /**
   * Find the minimum delta value that will give a non-zero score.
   * @return the minimum delta (likely to be &gt; 0 but this isn't guaranteed).
   */
  int minDelta();

  /**
   * Perform any end of life processing (including logging statistics)
   */
  void close();
}
