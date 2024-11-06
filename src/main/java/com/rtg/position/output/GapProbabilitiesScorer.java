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
import com.rtg.util.StringUtils;
import com.rtg.util.format.FormatInteger;
import com.rtg.util.format.FormatReal;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;



/**
 * Provides a distribution of scores for use when doing gapped merging.
 * Assumes that the i parameters are always multiples of step size.
 */
public class GapProbabilitiesScorer extends IntegralAbstract implements GapScorer {

  protected static double logProb(final double d) {
    assert d >= 0.0 && d <= 1.0 && !Double.isNaN(d);
    final double log = Math.log(d);
    assert log <= 0.0 && !Double.isNaN(log);
    return log;
  }

  protected final double[][] mScoreDistribution;

  protected final double[] mScoreMaxDistribution;

  protected final int mRowLength;

  /**
   * @param gapThreshold threshold
   * @param distribution the gap probability distribution for each offset
   * @param maxDistr maximum probability for each build offset.
   */
  GapProbabilitiesScorer(final double gapThreshold, final double[][] distribution, final double[] maxDistr) {
    mScoreDistribution = new double[distribution.length][];
    for (int i = 0; i < distribution.length; ++i) {
      final int length = distribution[i].length;
      mScoreDistribution[i] = new double[length];
      for (int j = 0; j < length; ++j) {
        mScoreDistribution[i][j] = logProb(distribution[i][j]);
      }
    }

    mScoreMaxDistribution = new double[mScoreDistribution.length];
    for (int i = 0; i < mScoreDistribution.length; ++i) {
      mScoreMaxDistribution[i] = logProb(maxDistr[i]);
    }

    mRowLength = distribution[0].length;
  }


  /**
   * Natural log of probability that the gap specified by the parameters will occur by chance.
   */
  @Override
  public double score(int buildStart, int buildEnd, int queryStart, int queryEnd) {
    final int buildLength = buildEnd - buildStart;
    final int queryLength = queryEnd - queryStart;
    //System.err.println("i=" + i + " j=" + j);
    if (buildLength < 0) {
      //overlap
      return buildLength == queryLength ? 0.0 : Double.NEGATIVE_INFINITY;
    }
    if (queryLength < 0) {
      return Double.NEGATIVE_INFINITY;
    }

    assert buildLength >= 0; // : "i=" + i + " io=" + io + " offset=" + mOffset + " stepSize=" + mStepSize;
    if (buildLength >= mScoreDistribution.length) {
      return Double.NEGATIVE_INFINITY;
    }
    final double[] d = mScoreDistribution[buildLength];
    if (queryLength >= d.length) {
      return Double.NEGATIVE_INFINITY;
    }
    final double log = d[queryLength];
    if (Double.isNaN(log)) {
      throw new IllegalStateException();
    }
    assert log <= 0.0 && !Double.isNaN(log);
    return log;
  }

  /**
   * Natural log of maximum probability that the gap specified by i (on the build) will occur by chance.
   */
  @Override
  public double scoreMax(int buildStart, int buildEnd, int queryStart, int queryEnd) {
    final int i = buildEnd - buildStart;
    if (i >= mScoreMaxDistribution.length) {
      return Double.NEGATIVE_INFINITY;
    }
    return mScoreMaxDistribution[i];
  }

  /**
   * Find the largest delta value that will give a non-zero probability.
   * @return the largest delta (likely to be &gt; 0 but this isn't guaranteed).
   */
  @Override
  public int maxDelta() {
    int m = Integer.MIN_VALUE;
    for (int is = 0; is < mScoreDistribution.length; ++is) {
      final int i = is;
      for (int j = 0; j < mRowLength; ++j) {
        final int delta = i - j;
        if (mScoreDistribution[is][j] > Double.NEGATIVE_INFINITY && delta > m) {
          m = delta;
        }
      }
    }
    return m;
  }

  /**
   * Find the smallest delta value that will give a non-zero probability.
   * @return the smallest delta (likely to be &lt; 0 but this isn't guaranteed).
   */
  @Override
  public int minDelta() {
    int m = Integer.MAX_VALUE;
    for (int is = 0; is < mScoreDistribution.length; ++is) {
      final int i = is;
      for (int j = mRowLength - 1; j >= 0; --j) {
        final int delta = i - j;
        if (mScoreDistribution[is][j] > Double.NEGATIVE_INFINITY && delta < m) {
          m = delta;
        }
      }
    }
    return m;
  }

  static final FormatInteger INDEX_FORMAT = new FormatInteger(3);
  static final FormatReal PROB_FORMAT = new FormatReal(2, 4);
  static final FormatReal SCORE_FORMAT = new FormatReal(3, 3);


  @Override
  public void toString(final StringBuilder sb) {
    sb.append(INDEX_FORMAT.blanks()).append("  ");
    for (int is = 0; is < mScoreDistribution.length; ++is) {
      final int i = is;
      sb.append("  [").append(INDEX_FORMAT.format(i)).append("]");
    }
    sb.append(StringUtils.LS);
    for (int j = 0; j < mRowLength; ++j) {
      sb.append("[").append(INDEX_FORMAT.format(j)).append("]");
      for (double[] aMScoreDistribution : mScoreDistribution) {
        final double sc = aMScoreDistribution[j];
        sb.append(sc == Double.NEGATIVE_INFINITY ? " -     " : SCORE_FORMAT.format(sc));
      }
      sb.append(StringUtils.LS);
    }
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (final double[] scoreDist : mScoreDistribution) {
      Exam.assertEquals(mRowLength, scoreDist.length);
      for (final double sc : scoreDist) {
        Exam.assertTrue(sc <= 0.0 && !Double.isNaN(sc));
      }
    }
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mScoreDistribution.length > 0);
    return true;
  }

  @Override
  public void setInternalBuildId(int internalId) {
  }

  @Override
  public void setQueryId(int internalId, Frame frame) {
  }


  @Override
  public void close() {
  }

}
