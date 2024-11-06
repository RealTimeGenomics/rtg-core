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
