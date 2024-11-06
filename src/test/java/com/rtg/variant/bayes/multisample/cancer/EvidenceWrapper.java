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
