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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;

import com.rtg.util.FlexArray;
import com.rtg.util.MathUtils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Contains transition table for homopolymer repeats as well as
 * limits of lengths where it will be applied. Variable names use:
 * <br>
 * R for read.
 * <br>
 * T for template.
 * <br>
 * The indices on the count an transition arrays are in the order:
 * <ul>
 * <li>nucleotide - range 0..3
 * <li>template length
 * <li>read length
 * </ul>
 */
public class HomoPolymerParams extends IntegralAbstract {

  /**
   * Read file containing counts for homopolymer transition table.
   * Each line has the form:
   * <br>
   *  tag i c_1 ... c_n
   *  <br>
   *  <ul>
   *  <li>tag is one of A+T or C+G
   *  <li>i the line number for the tag
   *  <li>counts starting at 1 up to some n.
   *  </ul>
   *
   * @param in the input file being read.
   * @return array of counts. The indices are: 0..1 - A+T or C+G, template length, read length.
   * @throws IOException whenever.
   */
  static int[][][] counts(final Reader in) throws IOException {
    final BufferedReader bin = new BufferedReader(in);
    @SuppressWarnings("unchecked")
    final FlexArray<int[]>[] temp = (FlexArray<int[]>[]) new FlexArray<?>[2];
    temp[0] = new FlexArray<>();
    temp[1] = new FlexArray<>();
    while (true) {
      //TODO make parsing and error handling more robust
      final String line = bin.readLine();
      if (line == null) {
        break;
      }
      final String[] split = line.split("\\t");
      if (split.length < 2) {
        throw new IOException("Invalid line in homopolymer calibration:" + line);
      }
      final int tag;
      if ("A+T".equals(split[0])) {
        tag = 0;
      } else if ("C+G".equals(split[0])) {
        tag = 1;
      } else {
        throw new IOException("Invalid line in homopolymer calibration:" + line);
      }
      final int t = Integer.parseInt(split[1]);
      final int[] cnt = new int[split.length - 2];
      for (int r = 0; r < cnt.length; ++r) {
        cnt[r] = Integer.parseInt(split[r + 2]);
      }
      temp[tag].set(t, cnt);
    }
    final int[][][] res = new int[2][][];
    res[0] = temp[0].toArray(new int[0][]);
    res[1] = temp[1].toArray(new int[0][]);
    return res;
  }

  /**
   * @param arithmetic use this for computing final possibilities.
   * @param counts raw counts used to compute probabilities.
   * @param complement iff true then compute 1.0 - probability rather than probability.
   * @return array of transition probabilities represented as possibilities.
   */
  static double[][][] transitions(final PossibilityArithmetic arithmetic, final int[][][] counts, final boolean complement) {
    final double[][][] res = new double[4][][];
    final double[][] at = transitions(arithmetic, counts[0], complement);
    final double[][] cg = transitions(arithmetic, counts[1], complement);
    //A and T share the same A+T statistics also C and G share the C+G statistics.
    res[0] = at;
    res[1] = cg;
    res[2] = cg;
    res[3] = at;
    return res;
  }

  private static double[][] transitions(final PossibilityArithmetic arithmetic, final int[][] counts, final boolean complement) {
    final double[][] res = new double[counts.length][];
    for (int i = 0; i < counts.length; ++i) {
      res[i] = transitions(arithmetic, counts[i], complement);
    }
    return res;
  }

  /**
   * Convert an array of counts to an array of probabilities represented as possibilities.
   * @param arithmetic use this for computing final possibilities.
   * @param counts raw counts used to compute probabilities.
   * @param complement iff true then compute 1.0 - probability rather than probability.
   * @return array of transition probabilities represented as possibilities.
   */
  static double[] transitions(final PossibilityArithmetic arithmetic, final int[] counts, final boolean complement) {
    if (counts == null) {
      return new double[0];
    }
    int total = 0;
    for (int count : counts) {
      total += count;
    }
    if (total == 0) {
      return new double[0];
    }
    final double[] res = new double[counts.length];
    for (int i = 0; i < counts.length; ++i) {
      final double tr = counts[i] / (double) total;
      final double t = complement ? 1.0 - tr : tr;
      res[i] = arithmetic.prob2Poss(t);
    }
    return res;
  }

  private final PossibilityArithmetic mArithmetic;

  private final int mMinR;

  private final int mMinT;

  private final double[][][] mTransitions;

  private final double[][][] mTransitionsC;

  /**
   * @param arithmetic possibility arithmetic - the transition results are returned in valid values for this arithmetic.
   * @param minR the shortest length of a repeat in a read which is treated as a homopolymer.
   * @param minT the shortest length of a repeat in a template which is treated as a homopolymer.
   * @param in contains text with arrays of counts for A+T and C+G
   * @throws IOException whenever.
   */
  public HomoPolymerParams(final PossibilityArithmetic arithmetic, final int minR, final int minT, final Reader in) throws IOException {
    this(arithmetic, minR, minT, counts(in));
  }

  //for private and testing use
  HomoPolymerParams(final PossibilityArithmetic arithmetic, final int minR, final int minT, final int[][][] counts) {
    this(arithmetic, minR, minT, transitions(arithmetic, counts, false), transitions(arithmetic, counts, true));
  }

  /**
   * @param arithmetic possibility arithmetic - the transition results are returned in valid values for this arithmetic.
   * @param minR the shortest length of a repeat in a read which is treated as a homopolymer.
   * @param minT the shortest length of a repeat in a template which is treated as a homopolymer.
   * @param transitions contains arrays of transition probabilities for A+T and C+G
   * @param complements as for transitions but containing the complement probabilities
   */
  HomoPolymerParams(final PossibilityArithmetic arithmetic, final int minR, final int minT, final double[][][] transitions,  final double[][][] complements) {
    mArithmetic = arithmetic;
    mMinR = minR;
    mMinT = minT;
    mTransitions = transitions;
    mTransitionsC = complements;
  }

  /**
   * @return the shortest length that is treated as a repeat on a read.
   */
  int minReadRepeat() {
    return mMinR;
  }

  /**
   * @return the shortest length that is treated as a repeat on the template.
   */
  int minTemplateRepeat() {
    return mMinT;
  }

  /**
   * Get a transition probability.
   * @param nt the current nucleotide (0=A ... 3=T)
   * @param t the length of the repeat on the template.
   * @param r  the length of the repeat on the read.
   * @return homopolymer transition probability as a possibility.
   */
  double transition(final int nt, final int t, final int r) {
    //System.err.println("nt=" + nt + " t=" + t + " r=" + r);
    if (r < mMinR || t < mMinT || t >= mTransitions[nt].length || r > mTransitions[nt][t].length) {
      return mArithmetic.zero();
    }
    return mTransitions[nt][t][r - 1];
  }

  /**
   * Get the complement of a transition probability.
   * @param nt the current nucleotide (0=A ... 3=T)
   * @param t the length of the repeat on the template.
   * @param r  the length of the repeat on the read.
   * @return (1.0 minus homopolymer transition probability) as a possibility.
   */
  double transitionC(final int nt, final int t, final int r) {
    //System.err.println("nt=" + nt + " t=" + t + " r=" + r);
    if (r < mMinR || t < mMinT || t >= mTransitionsC[nt].length || r > mTransitionsC[nt][t].length) {
      return mArithmetic.one();
    }
    return mTransitionsC[nt][t][r - 1];
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (int nt = 0; nt < 4; ++nt) {
      Exam.assertEquals(mTransitions[nt].length, mTransitionsC[nt].length);
      for (int t = 0; t < mTransitions[nt].length; ++t) {
        Exam.assertEquals(mTransitions[nt][t].length, mTransitionsC[nt][t].length);
        double total = 0.0;
        for (int r = 0; r < mTransitions[nt][t].length; ++r) {
          final double v = mTransitions[nt][t][r];
          total = mArithmetic.add(total, v);
          Exam.assertTrue(mArithmetic.isValidPoss(v));
          final double vc = mTransitionsC[nt][t][r];
          Exam.assertTrue(mArithmetic.isValidPoss(vc));
          Exam.assertTrue(MathUtils.approxEquals(mArithmetic.poss2Prob(v), 1.0 - mArithmetic.poss2Prob(vc), 0.000001));
        }
        if (mTransitions[nt][t].length > 0) {
          final double tp = mArithmetic.poss2Prob(total);
          Exam.assertTrue("tp=" + tp, MathUtils.approxEquals(tp, 1.0, 0.000001));
        }
      }
    }
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mArithmetic);
    Exam.assertTrue(mMinR >= 1);
    Exam.assertTrue(mMinT >= 1);
    Exam.assertEquals(4, mTransitions.length);
    Exam.assertEquals(4, mTransitionsC.length);
    return true;
  }
}
