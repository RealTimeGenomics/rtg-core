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
}
