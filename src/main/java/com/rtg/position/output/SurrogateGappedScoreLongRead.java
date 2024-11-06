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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.mode.Frame;
import com.rtg.reader.NamesInterface;

/**
 * Surrogate for <code>GappedScoreRegion</code>.
 *
 */
@TestClass(value = {"com.rtg.ngs.NgsLongTest"})
class SurrogateGappedScoreLongRead implements SurrogateRegion {


  private int mQueryId;
  private double mScore;
  private double mScoreThreshold;

  final SurrogateRegion initialize(final GappedScoreLongRead region, final int queryId, final Frame queryFrame, final double endCorr) {
    mQueryId = queryId;
    final double score = region.score();
    //System.err.println("score=" + score + "  endCorr=" + endCorr);
    mScore = score + endCorr;
    mScoreThreshold = region.mScoreThreshold;
    return this;
  }

  SurrogateGappedScoreLongRead(final GappedScoreLongRead region, final int queryId, final Frame[] subjectFrames, final Frame queryFrame, final double endCorr) {
    initialize(region, queryId, queryFrame, endCorr);
  }

  @Override
  public boolean write(final Appendable out, final NamesInterface subjectNames, final NamesInterface queryNames) {
    throw new UnsupportedOperationException();
  }

  @Override
  public boolean scoreAllowed() {
    return mScore >= mScoreThreshold;
  }

  /**
   * @see com.rtg.position.output.SurrogateRegion#writeHeader(java.lang.Appendable)
   */
  @Override
  public void writeHeader(final Appendable out) {
    throw new UnsupportedOperationException();
  }

  @Override
  public double score() {
    return -mScore;
  }

  @Override
  public int queryId() {
    return mQueryId;
  }
}
