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
    for (byte nt = 1; nt <= 4; ++nt) {
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
