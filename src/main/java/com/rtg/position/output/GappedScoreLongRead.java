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

import java.util.Locale;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.mode.Frame;
import com.rtg.util.array.ImmutableIntArray;
import com.rtg.util.integrity.Exam;

/**
 * Scores a read including correction for the ends of a read.
 * Also assumes that the residues in a gap hit (and uses the
 * probability correction in a gap).
 */
@TestClass(value = {"com.rtg.ngs.NgsLongTest"})
public class GappedScoreLongRead extends AbstractGappedRegion<GappedScoreLongRead> {

  protected static final double LN4 = Math.log(4.0);
  private final ImmutableIntArray mBuildLengths;
  protected double mScore = Double.NaN;
  protected final double mScoreIncrement;
  final double mScoreThreshold;

  /**
   * @param id an identifier unique within a particular <code>GappedOutput</code> object.
   * @param params parameter
   * @param distribution probability distribution for gap probabilities.
   * @param buildLengths build lengths
   */
  GappedScoreLongRead(final int id, final PositionParams params, final GapScorer distribution, final ImmutableIntArray buildLengths) {
    super(id, params.build(), distribution);
    mBuildLengths = buildLengths;
    final PositionOutputParams op = params.output();
    mScoreIncrement = LN4;
    final Double scoreThreshold = op.scoreThreshold();
    if (scoreThreshold == null) {
      mScoreThreshold = Double.NEGATIVE_INFINITY;
    } else {
      mScoreThreshold = scoreThreshold;
    }
  }

  @Override
  SurrogateRegion surrogateOutput(final SurrogateRegion surrogate, final int queryId, final Frame queryFrame, final int queryLength, final int queryEffectiveLength) {
    assert isInitialized();
    final int numberSFrames = mSubjectFrames.length;
    final int seqId = sequenceId() / numberSFrames;
    final int buildLength = mBuildLengths.get(seqId);
    final double endCorr = scoreEnd(buildLength, queryEnd(), queryLength);
    final SurrogateRegion res;
    if (surrogate == null) {
      res = new SurrogateGappedScoreLongRead(this, queryId, mSubjectFrames, queryFrame, endCorr);
    } else {
      res = ((SurrogateGappedScoreLongRead) surrogate).initialize(this, queryId, queryFrame, endCorr);
    }
    write();
    return res;
  }

  /**
   * Compute correction to score from ends of read.
   * @param buildLength length of sequence in build.
   * @param queryEnd last position that has been hit in the query.
   * @param queryLength total length of the query sequence.
   * @return correction to be added to total score.
   */
  double scoreEnd(final int buildLength, final int queryEnd, final int queryLength) {
    mDist.setInternalBuildId(mBuildSeqId);
    //System.err.println("scoreEnd buildLength=" + buildLength + " mBuildEnd=" + mBuildEnd + " queryLength=" + queryLength + " mQueryEnd=" + mQueryEnd);
    final double l1 = se(mDist, 0, mBuildStart, 0, mQueryStart);
    final double l2 = se(
        mDist,
        mBuildEnd + 1, buildLength,
        mQueryEnd + 1, queryLength
    );
    return l1 + l2;
  }

  /**
   * Compute scores for tails at start and end of read alignment.
   * @param dist probability distribution.
   * @param buildStart the start of the current build hit region
   * @param buildEnd the end of the previous build hit region
   * @param queryStart the start of the current query hit region
   * @param queryEnd the end of the previous query hit region
   * @return score for tail.
   */
  static double se(final GapScorer dist, final int buildStart, final int buildEnd, final int queryStart, final int queryEnd) {
    final double r;
    r = dist.scoreMax(buildStart, buildEnd, queryStart, queryEnd);
    return Double.isNaN(r) ? Double.NEGATIVE_INFINITY : r;
  }

  double score() {
    return mScore;
  }

  String formattedScore() {
    assert !Double.isNaN(score());
    return String.format(Locale.ROOT, "%1$04.2f", score());
  }

  @Override
  public void toString(final StringBuilder sb) {
    if (!isValid()) {
      sb.append("empty");
      return;
    }
    sb.append(mId);
    sb.append("[").append(mBuildStart).append("..").append(mBuildEnd).append("]");
    sb.append(":");
    sb.append("[").append(mQueryStart).append("..").append(queryEnd()).append("]");
    sb.append(formattedScore());
    sb.append("^").append(bucket());
    sb.append("->").append(next().mId);
  }

  @Override
  void initialize(final int seqId, final int buildStart, final int queryStart, final GapBuckets<GappedScoreLongRead> gb, final GappedScoreLongRead next, final boolean reverseFrame) {
    super.initialize(seqId, buildStart, queryStart, gb, next, reverseFrame);
    mScore = 0.0;
  }

  @Override
  void reset() {
    super.reset();
    mScore = Double.NaN;
  }

  @Override
  void merge(final GappedScoreLongRead that) {
    final double sc = gapScore(that);
    assert !Double.isNaN(sc) && sc <= 0.0;
    //System.err.println("merge this=" + this + " that=" + that + " GappedScore");
    if (this.mBuildStart > that.mBuildEnd) {
      //no overlap
      //System.err.println("increment sc=" + sc + " GappedScore");
      mScore += sc;
    }
    mScore += that.mScore;
    super.merge(that);
    //System.err.println("merge this=" + this + " mScoreIncrement=" + mScoreIncrement + " buildLength=" + buildLength() + " GappedScore");
  }

  @Override
  public boolean integrity() {
    super.integrity();
    if (isInitialized()) {
      Exam.assertTrue(mScore <= 0.0 && Double.isFinite(mScore));
    } else {
      Exam.assertTrue(Double.isNaN(mScore));
    }
    return true;
  }

  static final GappedRegionFactory<GappedScoreLongRead> FACTORY = new GappedRegionFactory<GappedScoreLongRead>() {
    @Override
    public GappedScoreLongRead region(int id, GapScorer distribution, PositionParams params, ImmutableIntArray buildLengths) {
      return new GappedScoreLongRead(id, params, distribution, buildLengths);
    }
  };
}
