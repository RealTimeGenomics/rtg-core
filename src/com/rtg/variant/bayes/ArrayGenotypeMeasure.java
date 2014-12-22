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
 * Measure backed by an array.
 */
public class ArrayGenotypeMeasure extends AbstractGenotypeMeasure {
  private final PossibilityArithmetic mArithmetic;
  private final double[] mMeasures;

  /**
   * Construct a measure backed by an array
   * @param arithmetic the arithmetic space
   * @param measures the values in the measure
   * @param hypotheses the hypotheses this measure covers
   */
  public ArrayGenotypeMeasure(PossibilityArithmetic arithmetic, double[] measures, Hypotheses<?> hypotheses) {
    super(hypotheses);
    mMeasures = measures;
    mArithmetic = arithmetic;
  }

  @Override
  public PossibilityArithmetic arithmetic() {
    return mArithmetic;
  }

  @Override
  public double measure(int hypothesis) {
    assert hypothesis < mMeasures.length;
    return mMeasures[hypothesis];
  }

  @Override
  public int size() {
    return mMeasures.length;
  }
}
