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
import com.rtg.visualization.AnsiDisplayHelper;
import com.rtg.visualization.DisplayHelper;

/**
 * Calculates a banded matrix of alignment probabilities.
 * This is similar to the <code>Gotoh</code> algorithm, but calculates the
 * probability of all paths to a given point, rather than just the best path.
 *
 * This class and all its subclasses are capable of using several different
 * representations of probabilities, for example, native doubles (0.0 to 1.0)
 * or natural logs.  We use the term 'Possibility' for one of these encoded
 * probability values, so all variables whose name ends with 'Poss' should
 * only be manipulated via one of the possibility methods.
 *
 */
@TestClass("com.rtg.variant.realign.ScoreMatrixTest")
public abstract class AbstractAllPaths extends IntegralAbstract implements AllPaths {

  private static final int FIELD_DP = 0;
  private static final int FIELD_WIDTH = 3 + FIELD_DP + (FIELD_DP > 0 ? 1 : 0);

  protected final RealignParams mParams;
  protected final PossibilityArithmetic mArith;

  // some constant/parameter probabilities.
  protected final double mZeroPoss;
  protected final double mOnePoss;
  protected final double mOneInFourPoss;
  protected final double mMatchPoss;
  protected final double mMisMatchPoss;
  protected final double mDeleteOpenPoss;
  protected final double mDeleteExtendPoss;
  protected final double mInsertOpenPoss;
  protected final double mInsertExtendPoss;
  protected final double mOneMinusDeleteExtendPoss;
  protected final double mOneMinusDeleteOpenPoss;
  protected final double mOneMinusInsertExtendPoss;
  protected final double mOneMinusDeleteInsertOpenPoss;
  protected final double mDeleteOpenInFourPoss;

  private final double mDelOpen;


  /**
   * The height of the matrix less 1 (so valid rows are from 0 .. <code>mLength</code> inclusive)
   */
  protected int mLength;

  /**
   * The minimum width of the matrix (i.e., the width of the first row)
   */
  protected int mWidth;

  protected int mMaxWidth;  //the current full width of the matrix
  protected int mMaxLength;  //the current full length of the matrix

  protected Environment mEnv;
  protected double mDeleteStartPoss;
  protected double mMatchStartPoss;

  private double[][] mMatch;  // Possibility values
  private double[][] mInsert; // Possibility values
  private double[][] mDelete; // Possibility values

  /**
   * Cumulative sums along the 'final' row (which is the first row for a reverse matrix).
   * These sums are the Possibilities that the read ends after
   * the given position.
   */
  protected double[] mEndScores;

  /**
   * A score matrix with the given maximum band width.
   *
   * @param arith helper object that does the arithmetic so that this code can be independent of the representation.
   * @param params the machine error model and related parameters.
   */
  protected AbstractAllPaths(final PossibilityArithmetic arith, final RealignParams params) {
    mArith = arith;
    //System.err.println(env.toString());
    mParams = params;
    mZeroPoss = mArith.zero();
    mOneInFourPoss = mArith.prob2Poss(0.25);
    mOnePoss = mArith.one();

    mMatchPoss = mArith.ln2Poss(mParams.matchLn());
    mMisMatchPoss = mArith.ln2Poss(mParams.misMatchLn());
    mDelOpen = Math.exp(mParams.deleteOpenLn());
    mDeleteOpenPoss = mArith.ln2Poss(mParams.deleteOpenLn());
    mDeleteExtendPoss = mArith.ln2Poss(mParams.deleteExtendLn());
    mInsertOpenPoss = mArith.ln2Poss(mParams.insertOpenLn());
    mInsertExtendPoss = mArith.ln2Poss(mParams.insertExtendLn());
    mOneMinusDeleteExtendPoss = mArith.prob2Poss(1.0 - Math.exp(mParams.deleteExtendLn()));
    mOneMinusDeleteOpenPoss = mArith.prob2Poss(1.0 - mDelOpen);
    mOneMinusInsertExtendPoss = mArith.prob2Poss(1.0 - Math.exp(mParams.insertExtendLn()));
    mOneMinusDeleteInsertOpenPoss = mArith.prob2Poss(1.0 - mDelOpen - Math.exp(mParams.insertOpenLn()));
    mDeleteOpenInFourPoss = mArith.multiply(mArith.ln2Poss(mParams.deleteOpenLn()), mOneInFourPoss);

    mLength = -1;
    mMaxLength = -1;
    mWidth = -1;
    mMaxWidth = -1;
    mDeleteStartPoss = Double.NaN;
    mMatchStartPoss = Double.NaN;
    mDeleteStartPoss = mArith.multiply(mArith.prob2Poss(mDelOpen), mArith.prob2Poss(1.0));
    mMatchStartPoss = mArith.prob2Poss(1.0 - mDelOpen);
  }

  @Override
  public PossibilityArithmetic arithmetic() {
    return mArith;
  }

  @Override
  public void setEnv(final Environment env) {
    final int width = 2 * env.maxShift() + 1;
    final int readLength = env.readLength();
    if (readLength > mMaxLength || width > mMaxWidth) {
      resizeMatrix(readLength, width);
    }
    mLength = readLength;
    mWidth = width;
    mEnv = env;
    calculateProbabilities();
  }

  void resizeMatrix(int length, int width) {
    mMaxLength = length;
    mMaxWidth = width;
    mMatch = new double[mMaxLength + 1][mMaxWidth];
    mInsert = new double[mMaxLength + 1][mMaxWidth];
    mDelete = new double[mMaxLength + 1][mMaxWidth];
    mEndScores = new double[mMaxWidth];
  }

  /**
   * How far a given row is shifted along the template, relative to
   * the expected start position of the read.  For CG, this accounts for
   * the CG gaps - it is used to keep the middle of the alignment band
   * in line with the most common offsets of the CG gaps.
   *
   * @param row one-based read position.
   * @return an offset along the template, from where the read is expected to start.
   * Smallest value will be for row 0.
   * Largest value will be for the last row.
   * But the values for intermediate rows may move backwards
   * and forwards between those bounds.
   */
  protected int rowOffset(final int row) {
    return row - mEnv.maxShift() - 1;
  }

  /**
   * @param row one-based read position
   * @return true if the current row is the start of a new read fragment (only used in CG reads)
   */
  protected boolean isFragmentStart(final int row) {
    return false;
  }

  /**
   * Calculate all the probabilities in the matrix and record the
   * total probability of the read.
   */
  protected abstract void calculateProbabilities();

  protected final void calculateInitialRow(final int initRow, final double delete, final double match) {
    for (int j = 0; j < mWidth; ++j) {
      mDelete[initRow][j] = delete;
      mMatch[initRow][j] = match;
      mInsert[initRow][j] = mZeroPoss;
    }
  }

  /**
   * Calculate the match/mismatch cost at a given cell of the matrix.
   *
   * @param i zero-based read position
   * @param j zero-based column position (<code>0 .. mWidth - 1</code>)
   * @return natural log of the probability
   */
  protected final double matchEq(final int i, final int j) {
    final byte te = mEnv.template(templateIndex(i, j));
    if (te == 0) {
      return mOneInFourPoss;
    }
    return matchEqTe(i, te);
  }

  /**
   * @param i read position in internal co-ordinates (0 based).
   * @param j template position in internal co-ordinates (0 based)
   * @return relative template position (suitable for interrogating environment).
   */
  protected int templateIndex(final int i, final int j) {
    return rowOffset(i + 1) + j;
  }

  /**
   * Calculate the match/mismatch probability at a given row and for a specified nucleotide.
   *
   * @param i zero-based read position
   * @param nt nucleotide (assumed to be on template, 0=N...4=T).
   * @return probability represented as a possibility.
   */
  final double matchEqTe(final int i, final byte nt) {
    final byte re = mEnv.read(i);
    //System.err.println("         i=" + i + " j=" + j + " te=" + te + " re=" + re + " sum=" + sum);
    if (re == 0) {
      return mOneInFourPoss;
    }

    final double q = mEnv.quality(i);
    final double q3 = q / 3.0;
    final double match;
    final double misMatch;
    if (nt == re) {
      match = 1.0 - q;
      misMatch = q3;
    } else {
      match = q3;
      misMatch = (1.0 - q3) / 3.0;
    }
    final double x = mArith.multiply(mMatchPoss, mArith.prob2Poss(match));
    final double y = mArith.multiply(mMisMatchPoss, mArith.prob2Poss(misMatch));
    return mArith.add(x, y);
  }

  @Override
  public abstract double totalScoreLn();

  final int length() {
    return mLength;
  }

  final int width() {
    return mWidth;
  }

  /**
   * Gets a visual representation of the alignment, in PPM image format
   * @return the PPM image string
   */
  public final String toPpm() {
    int fragTot = 0;
    for (int row = 0; row <= mLength; ++row) {
      if (isFragmentStart(row)) {
        ++fragTot;
      }
    }
    double best = rescale(mArith.poss2Ln(mInsert[0][0]), 0);
    for (int row = 0; row <= mLength; ++row) {
      for (int col = 0; col < mWidth; ++col) {
        best = Math.max(best, rescale(mArith.poss2Ln(mInsert[row][col]), 0));
        best = Math.max(best, rescale(mArith.poss2Ln(mMatch[row][col]), 0));
        best = Math.max(best, rescale(mArith.poss2Ln(mDelete[row][col]), 0));
      }
    }
    System.err.println("Best: " + best);

    final StringBuilder sb = new StringBuilder();
    sb.append("P3\n");
    sb.append(mLength + mWidth + 2).append(' ').append(mLength + fragTot + 3).append('\n');
    sb.append("100\n");

    final int rowStart = rowOffset(0);
    final int rowEnd = rowOffset(mLength) + mWidth;

    ppmTemplateRow(sb, rowStart, rowEnd);
    for (int row = 0; row <= mLength; ++row) {
      if (isFragmentStart(row)) { // insert a horizontal line where the CG gap/overlap happens.
        ppmTemplateRow(sb, rowStart, rowEnd);
      }
      sb.append(PPM_DNA_COLS[row == 0 ? 0 : mEnv.read(row - 1)]); // Left hand read base column
      // indent the row, so we line up with the template
      for (int i = 0; i < rowOffset(row) - rowStart; ++i) {
        ppmColor(sb, 0, 0, 0);
      }
      for (int col = 0; col < mWidth; ++col) {
        ppmColor(sb,
          rescale(mArith.poss2Ln(mInsert[row][col]), 0, best),
          rescale(mArith.poss2Ln(mMatch[row][col]), 0, best),
          rescale(mArith.poss2Ln(mDelete[row][col]), 0, best));
      }
      for (int i = rowOffset(row) - rowStart + mWidth; i < mLength + mWidth; ++i) {
        ppmColor(sb, 0, 0, 0);
      }
      sb.append(PPM_DNA_COLS[row == 0 ? 0 : mEnv.read(row - 1)]); // Right hand read base column
      sb.append(LS);
    }
    ppmTemplateRow(sb, rowStart, rowEnd);
    return sb.toString();
  }

  private static final String[] PPM_DNA_COLS = {
    "  50  50  50 ", // N
    "   0 100   0 ", // A
    "   0   0 100 ", // C
    "  50   0  50 ", // G
    " 100   0   0 ", // T
  };

  private void ppmTemplateRow(StringBuilder sb, int rowStart, int rowEnd) {
    sb.append(PPM_DNA_COLS[0]);
    for (int pos = rowStart; pos < rowEnd; ++pos) {
      sb.append(PPM_DNA_COLS[mEnv.template(pos)]);
    }
    sb.append(PPM_DNA_COLS[0]);
    sb.append(LS);
  }

  private static void ppmColor(StringBuilder sb, double ins, double match, double del) {
    sb.append(StringUtils.padLeft(Utils.realFormat(ins, FIELD_DP), FIELD_WIDTH + 1));
    sb.append(StringUtils.padLeft(Utils.realFormat(match, FIELD_DP), FIELD_WIDTH + 1));
    sb.append(StringUtils.padLeft(Utils.realFormat(del, FIELD_DP), FIELD_WIDTH + 1));
    sb.append(" ");
  }

  private static double rescale(final double x, int edge) {
    //be careful about outputting -0.0
    return Double.isInfinite(x) ? edge : x == 0.0 ? 0.0 : -x;
  }
  private static double rescale(final double x, int edge, double best) {
    //be careful about outputting -0.0
    //return Double.isInfinite(x) ? edge : x == 0.0 ? 0.0 : -100 * x / best;
    return Double.isInfinite(x) ? edge : x == 0.0 ? 0.0 : 100 - -100 * x / best;
  }

  @Override
  public final void toString(final StringBuilder sb) {
    if (mLength == -1) {
      sb.append("ScoreMatrix uninitialised").append(LS);
      return;
    }
    final char[] dnaChars = DNA.valueChars();
    final DisplayHelper dh = new AnsiDisplayHelper();
    final int rowStart = rowOffset(0);
    final int rowEnd = rowOffset(mLength) + mWidth;
    final double[] itmp = new double[mWidth];
    final double[] mtmp = new double[mWidth];
    final double[] dtmp = new double[mWidth];
    printTemplateRow(sb, rowStart, rowEnd, mEnv, "ScoreMatrix", FIELD_WIDTH);
    for (int row = 0; row <= mLength; ++row) {
      if (isFragmentStart(row)) {
        // show a horizontal line where the CG gap/overlap happens.
        sb.append("      ");
        for (int i = 0; i < (FIELD_WIDTH * 3 + 1) * (rowOffset(row) - rowStart + mWidth); ++i) {
          sb.append('-');
        }
        sb.append(LS);
      }
      sb.append(StringUtils.padBetween("[", 5, row + "]"));
      sb.append(row == 0 ? ' ' : dnaChars[mEnv.read(row - 1)]);
      // indent the row, so we line up with the template
      for (int i = 0; i < rowOffset(row) - rowStart; ++i) {
        sb.append(StringUtils.padLeft("|", 3 * FIELD_WIDTH + 1));
      }
      int bi = -1, bm = -1, bd = -1;
      for (int col = 0; col < mWidth; ++col) {
        bi = update(itmp, mArith.poss2Ln(mInsert[row][col]), bi, col);
        bm = update(mtmp, mArith.poss2Ln(mMatch[row][col]), bm, col);
        bd = update(dtmp, mArith.poss2Ln(mDelete[row][col]), bd, col);
      }
      for (int col = 0; col < mWidth; ++col) {
        sb.append(cell(dh, format(itmp[col]), col == bi, DisplayHelper.MAGENTA));
        sb.append(cell(dh, format(mtmp[col]), col == bm, DisplayHelper.GREEN));
        sb.append(cell(dh, format(dtmp[col]), col == bd, DisplayHelper.CYAN));
        sb.append("|");
      }
      sb.append(LS);
    }
    printTemplateRow(sb, rowStart, rowEnd, mEnv, "", FIELD_WIDTH);
  }

  private static int update(double[] arr, double val, int bestIn, int col) {
    int best = bestIn;
    arr[col] = val;
    if (best == -1 || arr[col] > arr[best]) {
      best = col;
    }
    return best;
  }

  private static String cell(DisplayHelper dh, String text, boolean mark, int col) {
    return mark ? dh.decorateForeground(text, col) : text;
  }

  static String format(final double x) {
    return format(x, FIELD_WIDTH);
  }

  private static String format(final double x, final int fw) {
    //be careful about outputting -0.0
    final String fst = Double.isInfinite(x) && x < 0.0 ? "" : x == 0.0 ? Utils.realFormat(0.0, FIELD_DP) : Utils.realFormat(-x, FIELD_DP);
    return StringUtils.padLeft(fst, fw);
  }

  static void printTemplateRow(final StringBuilder sb, final int rowStart, final int rowEnd, final Environment env, final String msg, int fw) {
    // print the header row, showing the template.
    final char[] dnaChars = DNA.valueChars();
    sb.append(StringUtils.padBetween(msg, 7 + 3 * fw, "|"));
    for (int pos = rowStart + 1; pos <= rowEnd; ++pos) {
      sb.append(dnaChars[env.template(pos)]);
      sb.append(StringUtils.padLeft(Integer.toString(env.absoluteTemplatePosition(pos)), 2 * fw - 1));
      sb.append(StringUtils.padLeft("|", fw + 1));
    }
    sb.append(LS);
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (int i = 0; i < mLength; ++i) {
      // the three matrices have exactly the same shape.
      Exam.assertEquals(mMaxWidth, mMatch[i].length);
      Exam.assertEquals(mMaxWidth, mDelete[i].length);
      Exam.assertEquals(mMaxWidth, mInsert[i].length);
      for (int j = 0; j < mWidth; ++j) {
        final double vi = mInsert[i][j];
        Exam.assertTrue("i=" + i + " j=" + j + " vi=" + vi, mArith.isValidPoss(vi));
        final double vm = mMatch[i][j];
        Exam.assertTrue("i=" + i + " j=" + j + " vm=" + vm, mArith.isValidPoss(vm));
        final double vd = mDelete[i][j];
        Exam.assertTrue("i=" + i + " j=" + j + " vd=" + vd, mArith.isValidPoss(vd));
      }
    }
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertEquals(mArith.zero(), mZeroPoss);
    Exam.assertEquals(mArith.one(), mOnePoss);
    Exam.assertTrue(mArith.isValidPoss(mOneInFourPoss));
    Exam.assertTrue(mArith.isValidPoss(mMatchPoss));
    Exam.assertTrue(mArith.isValidPoss(mMisMatchPoss));
    Exam.assertTrue(mArith.isValidPoss(mDeleteOpenPoss));
    Exam.assertTrue(mArith.isValidPoss(mDeleteExtendPoss));
    Exam.assertTrue(mArith.isValidPoss(mInsertOpenPoss));
    Exam.assertTrue(mArith.isValidPoss(mInsertExtendPoss));
    Exam.assertTrue(mArith.isValidPoss(mOneMinusDeleteExtendPoss));
    Exam.assertTrue(mArith.isValidPoss(mOneMinusDeleteOpenPoss));
    Exam.assertTrue(mArith.isValidPoss(mOneMinusInsertExtendPoss));
    Exam.assertTrue(mArith.isValidPoss(mOneMinusDeleteInsertOpenPoss));
    Exam.assertTrue(mArith.isValidPoss(mDeleteOpenInFourPoss));
    Exam.assertTrue(0.0 <= mDelOpen && mDelOpen <= 1.0 && !Double.isNaN(mDelOpen));

    if (mLength == -1) {
      Exam.assertEquals(-1, mWidth);
      Exam.assertTrue(mEnv == null);
      Exam.assertTrue(mMatch == null);
      Exam.assertTrue(mDelete == null);
      Exam.assertTrue(mInsert == null);
      Exam.assertTrue(Double.isNaN(mDeleteStartPoss));
      Exam.assertTrue(Double.isNaN(mMatchStartPoss));
    } else {
      Exam.assertTrue(mArith.isValidPoss(mDeleteStartPoss));
      Exam.assertTrue(mArith.isValidPoss(mMatchStartPoss));
      Exam.assertTrue(0 < mLength);
      Exam.assertTrue(0 < mWidth);
      Exam.assertTrue(mWidth <= mMaxWidth);
      Exam.assertTrue(mLength <= mMaxLength);
      Exam.assertTrue(rowOffset(0) < rowOffset(mLength));
      Exam.assertEquals(mMaxLength + 1, mMatch.length);
      Exam.assertEquals(mMaxLength + 1, mDelete.length);
      Exam.assertEquals(mMaxLength + 1, mInsert.length);
      Exam.assertNotNull(mEnv);
    }
    return true;
  }

  final double insert(final int row, final int col) {
    return mInsert[row][col];
  }

  final void setInsert(final int row, final int col, final double poss) {
//    assert !Double.isNaN(poss) && poss != Double.POSITIVE_INFINITY : poss + " @ " + row + ":" + col;
    mInsert[row][col] = poss;
  }

  final double delete(final int row, final int col) {
    return mDelete[row][col];
  }

  final void setDelete(final int row, final int col, final double poss) {
//    assert !Double.isNaN(poss) && poss != Double.POSITIVE_INFINITY : poss + " @ " + row + ":" + col;
    mDelete[row][col] = poss;
  }

  final double match(final int row, final int col) {
    return mMatch[row][col];
  }

  final void setMatch(final int row, final int col, final double poss) {
//    assert !Double.isNaN(poss) && poss != Double.POSITIVE_INFINITY : poss + " @ " + row + ":" + col;
    mMatch[row][col] = poss;
  }
}
