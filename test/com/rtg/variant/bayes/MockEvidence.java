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

package com.rtg.variant.bayes;

/**
 */
public class MockEvidence extends Evidence {
  private final double[] mProbs;
  private final int mRead;

  /**
   * @param description used for hypotheses etc.
   * @param notMap probability that the read doesn't map to this position.
   * @param prob probabilities for various nt.
   * @param read the hypothesis corresponding to the sequence (often a single nucleotide)
   * being used to update the model.
   */
  public MockEvidence(Description description, double notMap, final double[] prob, final int read) {
    super(description, notMap);
    mProbs = prob;
    mRead = read;
  }

  @Override
  public double probability(int index) {
    return mProbs[index];
  }

  @Override
  public double pe() {
    return 0.25;
  }

  @Override
  public int read() {
    return mRead;
  }

  @Override
  public double error() {
    return 0.1;
  }

  @Override
  public int getReadBasesLeft() {
    return 0;
  }

  @Override
  public int getReadBasesRight() {
    return 0;
  }

  @Override
  public void setReadBasesLeft(int readBasesLeft) {
    throw new UnsupportedOperationException();
  }

  @Override
  public void setReadBasesRight(int readBaseRight) {
    throw new UnsupportedOperationException();
  }

  @Override
  public boolean isForward() {
    return true;
  }

  @Override
  public boolean isReadPaired() {
    return false;
  }

  @Override
  public boolean isMated() {
    return false;
  }

  @Override
  public boolean isUnmapped() {
    return false;
  }

}
