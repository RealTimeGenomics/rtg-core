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

package com.rtg.variant.realign;

import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

/**
 */
public class ScoreFastUnderflowHomopolymer extends AbstractAllPathsFastUnderflow {

  private final HomoPolymerParams mHomopolymerParamsSimple;

  private final HomoPolymerParams mHomopolymerParamsLog;


  /**
   * @param params the machine error model and related parameters.
   * @param homoparams transition probabilities for homopolymer repeats (using simple arithmetic).
   * @param homoparamsLog  transition probabilities for homopolymer repeats (using log arithmetic).
   */
  public ScoreFastUnderflowHomopolymer(final RealignParams params, final HomoPolymerParams homoparams, HomoPolymerParams homoparamsLog) {
    super(params);
    mHomopolymerParamsSimple = homoparams;
    mHomopolymerParamsLog = homoparamsLog;
  }

  @Override
  protected AllPaths makeMatrix(PossibilityArithmetic arith, RealignParams params) {
    if (arith == SimplePossibility.SINGLETON) {
      return new HomopolymerMatrix(arith, params, mHomopolymerParamsSimple);
    } else if (arith == LogApproximatePossibility.SINGLETON) {
      return new HomopolymerMatrix(arith, params, mHomopolymerParamsLog);
    } else {
      throw new RuntimeException(arith.toString());
    }
  }
}
