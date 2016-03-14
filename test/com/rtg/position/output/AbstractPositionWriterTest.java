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
    final ISequenceParams seqParams = new MockSequenceParams(new MockReaderParams(new MockSequencesReader(SequenceType.DNA), SequenceMode.UNIDIRECTIONAL));
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
