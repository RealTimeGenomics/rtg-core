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
package com.rtg.metagenomics;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.metagenomics.matrix.Vector;

/**
 * Container for results from a species invocation
 */
@TestClass("com.rtg.metagenomics.SubBlockResultTest")
public abstract class BlockResult {
  protected final Vector mLikelihoods;
  private final Vector mR;
  private final Vector mVarianceLog;

  /**
   * @param r vector containing the species results for this block
   * @param varianceLog vector containing the variance for the species result for this block
   * @param likelihoods vector in which species likelihoods will be stored.
   */
  public BlockResult(Vector r, Vector varianceLog, Vector likelihoods) {
    mR = r;
    mVarianceLog = varianceLog;
    mLikelihoods = likelihoods;
  }

  /**
   * @return vector containing the species results for this block
   */
  public Vector getR() {
    return mR;
  }

  /**
   * @return vector containing the variance for the species result for this block
   */
  public Vector getVarianceLog() {
    return mVarianceLog;
  }

  /**
   * Set a result value
   * @param i the index of the result
   * @param r result value
   */
  public void setR(int i, double r) {
    mR.set(i, r);
  }

  /**
   * Set a result value
   * @param i the index of the result
   * @param varianceLog log of result variance
   */
  public void setVariance(int i, double varianceLog) {
    mVarianceLog.set(i, varianceLog);
  }

  /**
   * @return vector containing the likelihoods for the species result for this block
   */
  public Vector getLikelihoods() {
    return mLikelihoods;
  }

  /**
   * Set a result value
   * @param i the index of the result
   * @param likelihood the result likelihood
   */
  public void setLikelihood(int i, double likelihood) {
    mLikelihoods.set(i, likelihood);
  }
}
