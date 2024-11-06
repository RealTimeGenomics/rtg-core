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

import java.util.Arrays;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Calculates the probability distribution for a given position on the template,
 * by adding the probabilities of all possible alignments of the read. Uses a
 * banded <code>GOTOH</code>-like alignment algorithm in the forward and reverse
 * directions, then just sums the probabilities down a column.
 *
 */
@TestClass({"com.rtg.variant.realign.DeltaImplementationTest", "com.rtg.variant.realign.DeltaImplementationSimpleTest"})
public class DeltaImplementation extends Delta {

  protected final double mLog25;

  /**
   * Create a re-aligner that sums all probabilities of all alignments.
   * @param arith use for all internal arithmetic.
   * @param params parameters
   */
  public DeltaImplementation(final PossibilityArithmetic arith, final RealignParams params) {
    super(arith, new ScoreMatrix(arith, params), new ScoreMatrixReverse(arith, params));
    mLog25 = Math.log(0.25);
  }

  /**
   * Note that the result array is reused on each call, so clone
   * it if you must keep it.
   */
  @Override
  public double[] probabilitiesLn(final int templatePosition) {
    final int index = templatePosition - mEnv.absoluteTemplatePosition(0) + mEnv.maxShift();
    final int length = mEnv.readLength();
    final int width = mForward.width();
    if (index < 0 || index >= length + width - 1) {
      //before start or after end
      Arrays.fill(mResult, mLog25);
      return mResult;
    }
    //sum the parts not affected by the nucleotide changes
    double sum;
    if (index < width - 1) {
      //left
      //delete in the top column
      sum = mult(mForward.delete(0, index + 1), mReverse.delete(0, index + 1));
    } else {
      //all other cases
      sum = arithmetic().zero();
    }

    //vertical inserts
    final int readStart = Math.max(0, index - width + 2);
    final int readEnd = Math.min(length, index + 1);
    for (int readPos = readStart; readPos <= readEnd; ++readPos) {
      final double fwd = mForward.insert(readPos, index + 1 - readPos);
      final double rev = mReverse.insert(readPos, index + 1 - readPos);
      final double prod = mult(fwd, rev);
      sum = add(sum, prod);
    }

    //horizontal at ends
    //left
    for (int templatePos = index + 2; templatePos < width; ++templatePos) {
      final double fwdDel = mForward.delete(0, templatePos);
      final double revDel = mReverse.delete(0, templatePos);
      final double prodDel = mult(fwdDel, revDel);
      sum = add(sum, prodDel);

      final double fwdMa = mForward.match(0, templatePos);
      final double revMa = mReverse.match(0, templatePos);
      final double prodMa = mult(fwdMa, revMa);
      sum = add(sum, prodMa);
    }
    //right
    for (int templatePos = length; templatePos <= index; ++templatePos) {
      final double fwdDel = mForward.delete(length, templatePos - length);
      final double revDel = mReverse.delete(length, templatePos - length);
      final double prodDel = mult(fwdDel, revDel);
      sum = add(sum, prodDel);

      final double fwdMa = mForward.match(length, templatePos - length);
      final double revMa = mReverse.match(length, templatePos - length);
      final double prodMa = mult(fwdMa, revMa);
      sum = add(sum, prodMa);
    }

    Arrays.fill(mResult, sum);

    //System.err.println(" templatePosition=" + templatePosition + " index=" + index + " sum=" + sum);
    //sum the parts affected by the nucleotide changes
    for (int readPos = readStart; readPos <= readEnd; ++readPos) {
      final double rev = mReverse.match(readPos, index + 1 - readPos);
      if (readPos == 0) {
        final double fwd = mForward.match(0, index + 1);
        final double prod = mult(fwd, rev);
        //System.err.println(" templatePosition=" + templatePosition + " readPos=" + readPos + " fwd=" + fwd + " rev=" + rev);
        for (byte nt = 1; nt <= 4; ++nt) {
          //System.err.println(" nt=" + nt + " prod=" + prod);
          mResult[nt - 1] = add(mResult[nt - 1], prod);
        }
      } else {
        final double fwd = mForward.calculateMatch(readPos - 1, index - (readPos - 1));
        //System.err.println(" templatePosition=" + templatePosition + " readPos=" + readPos + " fwd=" + fwd + " rev=" + rev);
        for (byte nt = 1; nt <= 4; ++nt) {
          final double matchEq = mForward.matchEqTe(readPos - 1, nt);
          final double prod = mult(fwd, mult(matchEq, rev));
          //System.err.println(" nt=" + nt + " matchEq=" + matchEq + " prod=" + prod);
          mResult[nt - 1] = add(mResult[nt - 1], prod);
        }
      }
    }
    //    for (int nt = 1; nt <= 4; ++nt) {
    //      System.err.println(" templatePosition=" + templatePosition + " result[" + (nt - 1) + "]=" + mResult[nt - 1]);
    //    }
    for (int i = 0; i < mResult.length; ++i) {
      mResult[i] = arithmetic().poss2Ln(mResult[i]);
    }
    return mResult;
  }

}
