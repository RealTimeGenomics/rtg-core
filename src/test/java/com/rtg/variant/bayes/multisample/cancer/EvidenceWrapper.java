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
package com.rtg.variant.bayes.multisample.cancer;

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Evidence;

/**
 * Allow manually setting various evidence values
 */
public class EvidenceWrapper extends Evidence {
  private final double mQuality;
  private final double[] mProbabilities;
  private final int mRead;
  private final double mPE;

  /**
   * @param description description used by evidence
   * @param read value for read method
   * @param probabilities probabilities for probability method
   * @param quality value for error method
   * @param mapQ map error value
   * @param pe value for corresponding method
   */
  public EvidenceWrapper(Description description, int read, double[] probabilities, double quality, double mapQ, double pe) {
    super(description, mapQ);
    mProbabilities = probabilities;
    mQuality = quality;
    mRead = read;
    mPE = pe;
  }

  @Override
  public double probability(int index) {
    return mProbabilities[index];
  }

  @Override
  public double pe() {
    return mPE;
  }

  @Override
  public double error() {
    return mQuality;
  }

  @Override
  public int read() {
    return mRead;
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
  }

  @Override
  public void setReadBasesRight(int readBaseRight) {
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
  public boolean isFirst() {
    return true;
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
