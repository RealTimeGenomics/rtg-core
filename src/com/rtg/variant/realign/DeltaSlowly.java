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

import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * An inefficient realignment implementation that recalculates
 * the entire matrix for each nucleotide.
 *
 */
public class DeltaSlowly extends Delta {

  /** Forward score matrix used for temporary SNP calculations. */
  final AllPaths mForwardTmp;

  /**
   * @param arith used for all internal arithmetic.
   * @param forward Forward score matrix
   * @param forward2 Another forward score matrix, which is used for SNP calculations.
   * @param reverse Reverse score matrix
   */
  public DeltaSlowly(final PossibilityArithmetic arith, final ScoreMatrix forward, final ScoreMatrix forward2, final ScoreMatrixReverse reverse) {
    super(arith, forward, reverse);
    mForwardTmp = forward2;
  }

  @Override
  public double[] probabilitiesLn(final int index) {
    final int start = mEnv.absoluteTemplatePosition(0);
    for (byte nt = 1; nt <= 4; nt++) {
      final Environment snp = new EnvironmentSNP(mEnv, index - start, nt);
      mForwardTmp.setEnv(snp);
      // uncomment this code for a given template position to see all 4 combined matrices.
      //      if (index == 13 || index == 14) {
      //        mReverse.setEnv(snp); // make it match mForwardTmp
      //        System.out.println("================ nt = " + nt + "  total=" + mForwardTmp.totalScoreLn() + " ================");
      //        System.out.println(combinedTerse(mForwardTmp, mReverse));
      //        mReverse.setEnv(mEnv); // set it back again!
      //      }
      mResult[nt - 1] = mForwardTmp.totalScoreLn();
    }
    return mResult;
  }
}
