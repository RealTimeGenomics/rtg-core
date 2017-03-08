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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.mode.DNA;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Calculates the probability distribution for a given position on the
 * template, by adding the probabilities of all possible alignments of the read.
 *
 */
@TestClass({"com.rtg.variant.realign.DeltaImplementationTest", "com.rtg.variant.realign.DeltaImplementationSimpleTest"})
public abstract class Delta extends IntegralAbstract implements DeltaInterface {

  /**
   * Calculates alignment probabilities starting from the beginning of the read.
   */
  protected final ScoreMatrix mForward;

  /** Calculates alignment probabilities starting from the end of the read. */
  protected final ScoreMatrixReverse mReverse;

  protected Environment mEnv;

  private final PossibilityArithmetic mArithmetic;

  protected final double[] mResult = new double[4];

  /**
   * @param arith possibility arithmetic used in all calculations.
   * @param forward score matrix.
   * @param reverse score matrix.
   */
  protected Delta(final PossibilityArithmetic arith, final ScoreMatrix forward, final ScoreMatrixReverse reverse) {
    mArithmetic = arith;
    mForward = forward;
    mReverse = reverse;
  }

  @Override
  public void setEnv(final Environment env) {
    mForward.setEnv(env);
    mReverse.setEnv(env);
    mEnv = env;
  }

  @Override
  public void setReplace(String replace) {
    throw new UnsupportedOperationException();
  }

  @Override
  public final double totalScoreLn() {
    return mForward.totalScoreLn();
  }

  @Override
  public final double totalScore() {
    return mForward.totalScore();
  }

  @Override
  public final double total() {
    return mForward.total();
  }

  @Override
  public boolean underflow() {
    return mArithmetic.underflow(total());
  }

  @Override
  public PossibilityArithmetic arithmetic() {
    return mArithmetic;
  }

  /**
   * Compute the natural logs of the probabilities for each of the 4 possible nucleotides
   * at the specified template position.
   * @param index position on the template (0 based and in absolute co-ordinates independent of the
   * start position).
   * @return the probabilities, (raw, not normalized to sum to 1.0).
   */
  public abstract double[] probabilitiesLn(final int index);

  /**
   * Calculates the likelihood that the read covers a given template position.
   *
   * @param index an absolute template position.
   * @return the natural log of the probability.
   */
  public final double withinReadLn(final int index) {
    return mReverse.readStartsBeforeLn(index) + mForward.readEndsAfterLn(index);
  }

  protected static String formatTerse(final double x) {
    assert 0.0 <= x && x <= 1.0;
    final String res = Utils.realFormat(x, 3);
    if (res.equals("0.000")) {
      return " ";
    }
    if (res.equals("1.000")) {
      return "*";
    }
    if (res.startsWith("0.00")) {
      return ".";
    }
    if (res.startsWith("0.0")) {
      return "+";
    }
    if (x >= 0.9) {
      return "9";
    }
    final String res1 = Utils.realFormat(x, 1);
    return res1.substring(2, 3);
  }

  @Override
  public void toString(final StringBuilder sb) {
    if (mEnv == null) {
      sb.append("RealignImplementation uninitialised").append(LS);
      return;
    }
    sb.append("RealignImplementation(maxShift=").append(mEnv.maxShift()).append(")").append(LS);
    mForward.toString(sb);
    mReverse.toString(sb);
  }

  // <code>setEnv</code> needs to be called before calling this.
  String combinedToString() {
    final StringBuilder sb = new StringBuilder();
    combinedToString(sb);
    return sb.toString();
  }

  private static final int FIELD_DP = 3;
  private static final int FIELD_WIDTH = 3 + FIELD_DP + (FIELD_DP > 0 ? 1 : 0);

  private void combinedToString(final StringBuilder sb) {
    final char[] dnaChars = DNA.valueChars();
    final int rowStart = mForward.rowOffset(0);
    final int length = mForward.length();
    final int width = mForward.width();
    final int rowEnd = mForward.rowOffset(length) + width;
    ScoreMatrix.printTemplateRow(sb, rowStart, rowEnd, mEnv, "Combined", FIELD_WIDTH);
    final double total = mForward.total();
    for (int row = 0; row <= length; ++row) {
      sb.append(StringUtils.padBetween("[", 5, row + "]"));
      sb.append(row == 0 ? ' ' : dnaChars[mEnv.read(row - 1)]);
      // indent the row, so we line up with the template
      for (int i = 0; i < mForward.rowOffset(row) - rowStart; ++i) {
        sb.append(StringUtils.padLeft("|", 3 * FIELD_WIDTH + 1));
      }
      for (int j = 0; j < width; ++j) {
        final double ins = mArithmetic.divide(mArithmetic.multiply(mForward.insert(row, j), mReverse.insert(row, j)), total);
        sb.append(format(mArithmetic.poss2Prob(ins)));
        final double mat = mArithmetic.divide(mArithmetic.multiply(mForward.match(row, j),  mReverse.match(row, j)),  total);
        sb.append(format(mArithmetic.poss2Prob(mat)));
        final double del = mArithmetic.divide(mArithmetic.multiply(mForward.delete(row, j), mReverse.delete(row, j)), total);
        sb.append(format(mArithmetic.poss2Prob(del)));
        sb.append("|");
      }
      sb.append(LS);
    }
    ScoreMatrix.printTemplateRow(sb, rowStart, rowEnd, mEnv, "", FIELD_WIDTH);
  }

  private static String format(final double x) {
    assert 0.0 <= x && x <= 1.0;
    final String res = "  " + Utils.realFormat(x, FIELD_DP);
    if (res.equals("  0.000")) {
      return StringUtils.padLeft("", FIELD_WIDTH);
    }
    return res;
  }

  String combinedTerse(final Environment env) {
    final StringBuilder sb = new StringBuilder();
    combinedTerse(mArithmetic, env, sb, mForward, mReverse);
    return sb.toString();
  }

  static void combinedTerse(final PossibilityArithmetic arith, final Environment env, final StringBuilder sb, final AbstractAllPaths fwd, final AbstractAllPaths rev) {
    //final Environment env = fwd.getEnv();
    final char[] dnaChars = DNA.valueChars();
    final int rowStart = fwd.rowOffset(0);
    final int length = fwd.length();
    final int width = fwd.width();
    final int rowEnd = fwd.rowOffset(length) + width;
    printTemplateRowTerse(sb, rowStart, rowEnd, env);
    final double total = fwd.total();
    for (int row = 0; row <= length; ++row) {
      sb.append(StringUtils.padBetween("[", 5, row + "]"));
      sb.append(row == 0 ? ' ' : dnaChars[env.read(row - 1)]);
      // indent the row, so we line up with the template
      for (int i = 0; i < fwd.rowOffset(row) - rowStart; ++i) {
        sb.append("   |");
      }
      for (int j = 0; j < width; ++j) {
        final double ins = arith.divide(arith.multiply(fwd.insert(row, j), rev.insert(row, j)), total);
        sb.append(formatTerse(arith.poss2Prob(ins)));
        final double mat = arith.divide(arith.multiply(fwd.match(row, j),  rev.match(row, j)),  total);
        sb.append(formatTerse(arith.poss2Prob(mat)));
        final double del = arith.divide(arith.multiply(fwd.delete(row, j), rev.delete(row, j)), total);
        sb.append(formatTerse(arith.poss2Prob(del)));
        sb.append("|");
      }
      sb.append(LS);
    }
    printTemplateRowTerse(sb, rowStart, rowEnd, env);
  }

  static void printTemplateRowTerse(final StringBuilder sb, final int rowStart, final int rowEnd, final Environment env) {
    // print the header row, showing the template.
    final char[] dnaChars = DNA.valueChars();
    sb.append(StringUtils.spaces(9));
    for (int pos = rowStart + 1; pos <= rowEnd; ++pos) {
      sb.append(dnaChars[env.template(pos)]);
      sb.append("   ");
    }
    sb.append(LS);
  }
  @Override
  public boolean integrity() {
    Exam.assertNotNull(mForward);
    Exam.assertNotNull(mReverse);
    Exam.assertEquals(mForward.length(), mReverse.length());
    Exam.assertEquals(mForward.width(), mReverse.width());
    return true;
  }

  protected double add(final double a, final double b) {
    return mArithmetic.add(a, b);
  }

  protected double mult(final double a, final double b) {
    return mArithmetic.multiply(a, b);
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    if (mEnv == null) {
      return true;
    }
    final double total = totalScoreLn();
    final double epsilon = Math.abs(total / 10000);
    Exam.assertEquals(total, mReverse.totalScoreLn(), epsilon);
    Exam.assertEquals(total, mForward.totalScoreLn(), epsilon);
    // test that various cuts sum to the total.
    //The following tests will almost certainly not work for CG
    //horizontal cuts
    final int length = mForward.length();
    final int width = mForward.width();
    for (int i = 0; i <= length; ++i) {
      double sum = mArithmetic.zero();
      for (int j = 0; j < width; ++j) {
        sum = add(sum, mult(mForward.delete(i, j), mReverse.delete(i, j)));
        sum = add(sum, mult(mForward.match(i, j), mReverse.match(i, j)));
      }
      Exam.assertEquals(total, mArithmetic.poss2Ln(sum), epsilon);
    }
    //vertical cuts on the left
    //j = 0 is special case covered elsewhere - ignore it now
    for (int j = 1; j < width; ++j) {
      double sum = mArithmetic.zero();
      for (int i = 0; i <= j; ++i) {
        final int k = j - i;
        sum = add(sum, mult(mForward.insert(i, k), mReverse.insert(i, k)));
        sum = add(sum, mult(mForward.match(i, k), mReverse.match(i, k)));
      }
      sum = add(sum, mult(mForward.delete(0, j), mReverse.delete(0, j)));
      for (int i = j + 1; i < width; ++i) {
        sum = add(sum, mult(mForward.delete(0, i), mReverse.delete(0, i)));
        sum = add(sum, mult(mForward.match(0, i), mReverse.match(0, i)));
      }
      Exam.assertEquals(total, mArithmetic.poss2Ln(sum), epsilon);
    }
    //vertical cut in the middle
    for (int j = width, l = 1; j < length; ++j, ++l) {
      double sum = mult(mForward.delete(l, j - l), mReverse.delete(l, j - l));
      for (int i = l; i < l + width; ++i) {
        final int k = j - i;
        //System.err.println("i=" + i + " l=" + l + " j=" + j + " k=" + k);
        sum = add(sum, mult(mForward.insert(i, k), mReverse.insert(i, k)));
        sum = add(sum, mult(mForward.match(i, k), mReverse.match(i, k)));
      }
      Exam.assertEquals(total, mArithmetic.poss2Ln(sum), epsilon);
    }
    //vertical cuts on the right
    for (int j = length; j < width + length - 1; ++j) {
      double sum = mArithmetic.zero();
      for (int i = j - width + 1; i < length; ++i) {
        final int k = j - i;
        sum = add(sum, mult(mForward.insert(i, k), mReverse.insert(i, k)));
        sum = add(sum, mult(mForward.match(i, k), mReverse.match(i, k)));
      }
      final int l = length - 1;
      for (int i = l; i < j; ++i) {
        final int k = i - l;
        sum = add(sum, mult(mForward.delete(l, k), mReverse.delete(l, k)));
        sum = add(sum, mult(mForward.match(l, k), mReverse.match(l, k)));
      }
      Exam.assertEquals(total, mArithmetic.poss2Ln(sum), epsilon);
    }

    //diagonal cuts
    for (int j = width - 1; j < length; ++j) {
      double sum = mArithmetic.zero();
      for (int i = j - width + 1, l = j; i <= l; ++i, --l) {
        final int k = l - i;
        sum = add(sum, mult(mForward.insert(i, k), mReverse.insert(i, k)));
        sum = add(sum, mult(mForward.match(i, k), mReverse.match(i, k)));
        sum = add(sum, mult(mForward.delete(i, k), mReverse.delete(i, k)));
      }
      for (int i = j - width + 2, l = j; i <= l; ++i, --l) {
        final int k = l - i;
        sum = add(sum, mult(mForward.match(i, k), mReverse.match(i, k)));
      }
      Exam.assertEquals(total, mArithmetic.poss2Ln(sum), epsilon);
    }
    return true;
  }

}
