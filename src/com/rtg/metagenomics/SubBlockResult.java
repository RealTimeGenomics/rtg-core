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

import com.rtg.metagenomics.matrix.Vector;

/**
 */
public class SubBlockResult extends BlockResult {

  private final double mL;

  /**
   * @param r vector containing the species results for this block (one per member of the block)
   * @param varianceLog vector containing the variance for the species result for this block (one per global id)
   * @param likelihoods vector containing the likelihoods (one per global id)
   * @param l the final objective function value after solving
   */
  public SubBlockResult(Vector r, Vector varianceLog, Vector likelihoods, double l) {
    super(r, varianceLog, likelihoods);
    mL = l;
  }

  /**
   * @return the final objective function value after solving
   */
  public double getL() {
    return mL;
  }
}
