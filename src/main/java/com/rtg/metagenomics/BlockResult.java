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
