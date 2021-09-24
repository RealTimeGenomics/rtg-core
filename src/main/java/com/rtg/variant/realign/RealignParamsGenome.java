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
