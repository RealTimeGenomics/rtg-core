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

import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Implementation of GenotypeMeasure that returns Nan for non identity
 * Used by disease caller which doesn't have a valid non identity score or strictly speaking an actual GQ
 */
public class NoNonIdentityMeasure implements GenotypeMeasure {
  private final GenotypeMeasure mInternal;

  /**
   *
   * @param internal the measure to encapsulate
   */
  public NoNonIdentityMeasure(GenotypeMeasure internal) {
    mInternal = internal;
  }

  @Override
  public PossibilityArithmetic arithmetic() {
    return mInternal.arithmetic();
  }

  @Override
  public double measure(int hypothesis) {
    return mInternal.measure(hypothesis);
  }

  @Override
  public int size() {
    return mInternal.size();
  }

  @Override
  public double bestPosterior() {
    return mInternal.bestPosterior();
  }

  @Override
  public double nonIdentityPosterior() {
    return Double.NaN;
  }

  @Override
  public int best() {
    return mInternal.best();
  }

  @Override
  public int reference() {
    return mInternal.reference();
  }

  @Override
  public Hypotheses<?> hypotheses() {
    return mInternal.hypotheses();
  }

  @Override
  public String toString() {
    return "NoNonIdentityMeasure{" + "mInternal=" + mInternal + '}';
  }
}
