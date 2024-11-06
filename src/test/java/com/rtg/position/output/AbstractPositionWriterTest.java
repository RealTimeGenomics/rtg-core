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

import com.rtg.launcher.BuildParams;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.mode.Frame;
import com.rtg.mode.SequenceMode;
import com.rtg.mode.SequenceType;
import com.rtg.reader.MockSequencesReader;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 */
public class AbstractPositionWriterTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  protected static class MockGappedRegion extends AbstractGappedRegion<MockGappedRegion> {

    private final double mScore;

    public MockGappedRegion(final int id, final BuildParams params, final GapScorer dist, final double score) {
      super(id, params, dist);
      mScore = score;
    }

    @Override
    SurrogateRegion surrogateOutput(final SurrogateRegion surrogate, final int queryId, final Frame queryFrame, final int queryLength, final int queryEffectiveLength) {
      final SurrogateRegion res;
      if (surrogate == null) {
        res = new MockSurrogateGappedRegion(this, queryId, queryFrame, mSubjectFrames, mScore);
      } else {
        res =  ((MockSurrogateGappedRegion) surrogate).initialize(this, queryId, queryFrame, mScore);
      }
      mState = State.WRITTEN;
      return res;

    }
  }
  private static class MockSurrogateGappedRegion extends SurrogateGappedRegion {

    private double mScore;

    MockSurrogateGappedRegion(final MockGappedRegion region, final int queryId, final Frame queryFrame, final Frame[] subjectFrames, final double score) {
      super(region, queryId, queryFrame, subjectFrames);
      mScore = score;
    }

    public SurrogateRegion initialize(final MockGappedRegion region, final int queryId, final Frame queryFrame, final double score) {
      final SurrogateRegion ret = super.initialize(region, queryId, queryFrame);
      ((MockSurrogateGappedRegion) ret).mScore = score;
      return ret;
    }

    @Override
    public double score() {
      return mScore;
    }
  }

  protected AbstractGappedRegion<?> getRegion(final int id, final int stepSize, final int wordSize, final double score) throws IOException {
    final GappedDistribution distr = new GappedDistribution(stepSize, wordSize, GappedDistribution.distrParams(3 * stepSize));
    final GapScorer prob = distr.probabilities();
    final ISequenceParams seqParams = new MockSequenceParams(new MockReaderParams(new MockSequencesReader(SequenceType.DNA)), SequenceMode.UNIDIRECTIONAL);
    final BuildParams params = BuildParams.builder().windowSize(wordSize).stepSize(stepSize).sequences(seqParams).create();
    final MockGappedRegion gr = new MockGappedRegion(id, params, prob, score);
    final BuildParams bparams = GapBucketsTest.makeParams(8, new int[] {30, 31, 48}, 0, SequenceMode.UNIDIRECTIONAL);
    final GapBucketsInfo bucketInfo = new GapBucketsInfo(bparams, 1, 16, 1000);
    final GapBuckets<MockGappedRegion> bu = new GapBuckets<>(bucketInfo);
    gr.initialize(2, 1, 32, bu, gr, false);
    return gr;
  }

  public void testDummy() {
    //do nothing - here to shut the warnigns up
  }
}
