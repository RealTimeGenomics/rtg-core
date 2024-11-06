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

import com.rtg.util.diagnostic.SpyCounter;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;


/**
 * Calculates a banded matrix of alignment probabilities.
 * This is similar to the <code>Gotoh</code> algorithm, but calculates the
 * probability of all paths to a given point, rather than just the best path.
 *
 * This does just the forward direction matrix.
 *
 * Allowance is made for homopolymer transition probabilities.
 *
 */
public class HomopolymerMatrix extends ScoreMatrix {

  private static final SpyCounter SPY = new SpyCounter("HomopolymerMatrix");

  private final HomoPolymerParams mHomoparams;

  /**
   * A score matrix with the given maximum band width.
   * @param arith helper object that does the arithmetic so that this code can be independent of the representation.
   * @param params the machine error model and related parameters.
   * @param homoparams transition probabilities for homopolymer repeats.
   */
  public HomopolymerMatrix(final PossibilityArithmetic arith, final RealignParams params, final HomoPolymerParams homoparams) {
    super(arith, params);
    mHomoparams = homoparams;
    SPY.increment();
  }

  /**
   * Calculate the forward diagonal probabilities, starting from cell <code>(i,j)</code>.
   * This does NOT include the match/mismatch probability of the destination cell.
   *
   * @param i zero-based read position
   * @param j column
   * @return Possibility
   */
  private double calculateMatchCarefully(final int i, final int j) {
    if (j < 0 || i < 0) {
      return mArith.zero();
    }
    return calculateMatch(i, j);
  }

  //  private String pr(final double x) {
  //    return Utils.realFormat(-mArith.poss2Ln(x), 3);
  //  }

  private EnvironmentHomopolymer env() {
    return (EnvironmentHomopolymer) mEnv;
  }

  @Override
  protected final void matchIt(final int i, final int j) {
    final int ti = templateIndex(i - 1, j);
    final double cm = calculateMatch(i - 1, j);
    final double cmeq = mArith.multiply(cm, matchEq(i - 1, j));
    final double newMatch;
    final boolean same = mEnv.template(ti) == mEnv.read(i - 1) && mEnv.template(ti) != 0;
    //System.err.println("i=" + i + " j=" + j + " ti=" + ti + " t=" + mEnv.template(ti) + " r=" + mEnv.read(i - 1) + " s=" + same + " ts=" + mEnv.templateStart(ti) + " rs=" + mEnv.readStart(i - 1) + " te=" + mEnv.templateEnd(ti) + " re=" + mEnv.readEnd(i - 1) + " mt=" + mHomoparams.minTemplateRepeat() + " mr=" + mHomoparams.minReadRepeat());
    if (same && env().templateStart(ti) >= mHomoparams.minTemplateRepeat() && env().readStart(i - 1) >= mHomoparams.minReadRepeat()) {
      //System.err.println("Bingo1" + " i=" + i + " j=" + j);
      //start of a repeat on both axes
      //System.err.println(mHomoparams);
      final double transitionC = mHomoparams.transitionC(env().template(ti) - 1, env().templateStart(ti), env().readStart(i - 1));
      newMatch =  mArith.multiply(cmeq, transitionC);
      //System.err.println("trC=" + pr(transitionC) + " cmeq=" + pr(cmeq) + " nm=" + pr(newMatch));

    } else if (same && env().templateEnd(ti) >= mHomoparams.minTemplateRepeat() && env().readEnd(i - 1) >= mHomoparams.minReadRepeat()) {
      //end of a repeat on both axes
      //System.err.println("Bingo2" + " i=" + i + " j=" + j);
      final double transition = mHomoparams.transition(env().template(ti) - 1, env().templateEnd(ti), env().readEnd(i - 1));
      //System.err.println(mArith.poss2Prob(transition));
      final int ip = i - 1 - env().readEnd(i - 1);
      final int jp = j - env().templateEnd(ti);
      final double cmShortCircuit = calculateMatchCarefully(ip, jp);
      final double trcm = mArith.multiply(transition, cmShortCircuit);
      newMatch = mArith.add(cmeq, trcm);
    } else {
      //non-shortcircuit
      newMatch = cmeq;
    }
    setMatch(i, j, newMatch);
  }
}
