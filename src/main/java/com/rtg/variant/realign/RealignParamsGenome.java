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

import com.rtg.util.machine.MachineType;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Suitable for use when using all paths on human genomic template.
 * See <code>com/rtg/variant/priors/prior/human.properties</code> for where the numbers were taken from.
 * TODO make this suitable for use with any priors file.
 */
public final class RealignParamsGenome implements RealignParams {

  /** Single instance for human genome. */
  public static final RealignParamsGenome SINGLETON = new RealignParamsGenome();

  private static final double MISMATCH = 0.00071 + 0.00053;

  private static final double MISMATCH_LN = Math.log(MISMATCH);

  private static final double MATCH_LN = Math.log(1.0 - MISMATCH);

  private static final double INDEL_OPEN = 0.00015 / 2.0;

  private static final double INDEL_OPEN_LN = Math.log(INDEL_OPEN);

  private static final double FITTED_INDEL_EXTEND = 0.5;

  private static final double INDEL_EXTEND_LN = Math.log(FITTED_INDEL_EXTEND);

  private RealignParamsGenome() { }

  @Override
  public double insertOpenLn() {
    return INDEL_OPEN_LN;
  }

  @Override
  public double insertExtendLn() {
    return INDEL_EXTEND_LN;
  }

  @Override
  public double deleteOpenLn() {
    return INDEL_OPEN_LN;
  }

  @Override
  public double deleteExtendLn() {
    return INDEL_EXTEND_LN;
  }

  @Override
  public double matchLn() {
    return MATCH_LN;
  }

  /**
   * @return the probability of a mismatch (really the probability of a single nucleotide mutation on the genome).
   */
  public double misMatch() {
    return MISMATCH;
  }

  @Override
  public double misMatchLn() {
    return MISMATCH_LN;
  }

  @Override
  public MachineType machineType() {
    return null;
  }

  @Override
  public int gapStart(int gap) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int gapEnd(int gap) {
    throw new UnsupportedOperationException();
  }

  @Override
  public double gapFreqLn(int gap, int width) {
    throw new UnsupportedOperationException();
  }

  @Override
  public double[][] gapDistributionPoss(PossibilityArithmetic arith) {
    throw new UnsupportedOperationException();
  }

}
