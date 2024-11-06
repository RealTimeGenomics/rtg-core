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
