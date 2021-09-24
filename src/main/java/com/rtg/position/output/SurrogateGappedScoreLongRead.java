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
