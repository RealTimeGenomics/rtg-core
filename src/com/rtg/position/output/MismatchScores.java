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

import com.rtg.mode.Frame;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Scoring logic for long reads.
 */
public class MismatchScores extends IntegralAbstract implements GapScorer {

  private final int mMaxGap;

  /** Gaps less or equal to this are permitted to have indels. */
  private final int mMaxIndel;
  protected final int mStepSize;
  private final int mWordSize;
  private final int mGapConstant;

  private long mScoreCalls = 0;
  private long mScoreMaxCalls = 0;

  /**
   * Constructor
   * @param maxGap maximum allowed gap.
   * @param maxIndel Maximum indel we are looking for
   * @param wordSize window size
   * @param stepSize step size
   */
  public MismatchScores(final int maxGap, final int maxIndel, final int wordSize, final int stepSize) {
    mMaxGap = maxGap;
    mMaxIndel = maxIndel;
    mWordSize = wordSize;
    mStepSize = stepSize;
    mGapConstant = (mWordSize - 1) / mStepSize + 1;
  }

  protected boolean allowed(final int i, final int j) {
    if (i <= mMaxGap && j <= mMaxGap) {
      final int delta = i - j;
      if (delta <= mMaxIndel && -mMaxIndel <= delta) {
        return true;
      }
    }
    return false;
  }

  @Override
  public double scoreMax(int buildStart, int buildEnd, int queryStart, int queryEnd) {
    ++mScoreMaxCalls;
    return score(buildStart, buildEnd, buildStart, buildEnd);
  }

  /**
   * Find the minimum number of substitutions that explains the gap. <p>
   *
   * General formula is: <code>ceiling(m/l)</code>, where <code>m</code> is the number of hashes that
   * didn't match over the entire gap and <code>l</code> is the number of words thats mismatch
   * is explained by a substitution in the first step of the gap. <p>
   *
   * <code>m = ceiling((gapSize + wordSize / stepSize)) - 1</code><br>
   * <code>l = ceiling(wordSize / stepSize)</code><p>
   *
   * Formula to turn <code>ceiling()</code> into <code>floor()</code>: <p>
   * <code>ceiling(a / b) = floor((a - 1) / b) + 1</code><p>
   *
   * The formula assumes a single substitution will explain <code>l</code> hashes and there
   * will need be another substitution to explain the next <code>l</code> hashes and so
   * on until at least <code>m</code> are explained.
   *
   * @param gap the size of the gap
   * @return number of substitutions that explains the gap.
   */
  int minSubs(final int gap) {
    assert gap >= 0;
    final int m = (gap + mWordSize - 1) / mStepSize;
    final int s = (m  - 1) / mGapConstant + 1;
    assert s >= 0;
    return s;
  }

  /**
   * Return the minimum possible mismatch score that explains the gap described
   * by <code>i</code> and <code>j</code>
   *
   * @param buildStart the start of the current build hit region
   * @param buildEnd the end of the previous build hit region
   * @param queryStart the start of the current query hit region
   * @param queryEnd the end of the previous query hit region
   * @return the negative of the mismatch score (since the caller thinks bigger is better)
   */
  @Override
  public double score(int buildStart, int buildEnd, int queryStart, int queryEnd) {
    ++mScoreCalls;
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

    if (buildLength == 0 && queryLength == 0) {
      return 0;
    }
    assert buildLength >= 0 && queryLength >= 0 : "Step Size: " + mStepSize + " i: " + buildLength + " j: " + queryLength;
    if (!allowed(buildLength, queryLength)) {
      return Double.NEGATIVE_INFINITY;
    }
    final int res;
    final int min = minSubs(buildLength);
    if (buildLength == queryLength) {
      res = min;
    } else {
      //min is the minimum number of window misses that must be "explained away"
      final int del;
      if (buildLength < queryLength) {
        del = queryLength - buildLength;
      } else {
        del = buildLength - queryLength;
      }

      assert del >= 1;
      if (del > min) {
        res = 2 * del;
      } else {
        res = del + min;
      }
    }
    assert res >= 0;
    //Warning we return the NEGATIVE of the mismatch score because GappedOutput thinks bigger is better.
    return -res;
  }

  @Override
  public int maxDelta() {
    return mMaxIndel;
  }

  @Override
  public int minDelta() {
    return -mMaxIndel;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("MismatchProbabilities");
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mMaxGap >= 1);
    Exam.assertTrue(mMaxIndel >= 0);
    Exam.assertTrue(mGapConstant >= 1);
    Exam.assertTrue(mWordSize >= mStepSize);
    Exam.assertTrue(mStepSize >= 1);
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
    Diagnostic.developerLog("GappedOutput score function calls: " + mScoreCalls);
    Diagnostic.developerLog("GappedOutput scoreMax function calls: " + mScoreMaxCalls);
  }

}
