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
package com.rtg.alignment;

import java.util.Arrays;
import java.util.Locale;

import com.rtg.mode.DNA;
import com.rtg.reader.CgUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.format.FormatReal;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.realign.Environment;
import com.rtg.variant.realign.EnvironmentImplementation;
import com.rtg.variant.realign.InvertedEnvironment;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.RealignParamsImplementation;

/**
 * This finds good alignments for Complete Genomics reads (length = 35 or 29),
 * taking into consideration the probabilities of varying overlap and gap sizes.
 *
 * It implements the Needleman-Wunsch global alignment algorithm with the improvements
 * described by Osamu Gotoh for handling affine gap open penalties.
 * See Journal of Molecular Biology, Volume 162, Issue 3, 15 December 1982, Pages 705-708.
 * This handles affine gap penalties (a gap open penalty plus an extend penalty),
 * but still guarantees to find the lowest cost alignment.
 *
 * It constructs three 2D arrays of probabilities using dynamic programming:
 * <ul>
 *   <li><code>insert</code> is the probability of an alignment that ends with an insertion.</li>
 *   <li><code>match</code> is the probability of an alignment that ends with a diagonal.</li>
 *   <li><code>delete</code> is the probability of an alignment that ends with a deletion.</li>
 * </ul>
 * For consistency, we always refer to the three arrays in the above (clockwise) order:
 * insert scores come from the left (<code>9 am</code>), match scores come from the upper left
 * diagonal (<code>10:30 am</code>) and delete scores come from above (<code>12 am</code>).  The
 * dump display also follows this order.
 *
 */
public class CgGotohEditDistance extends IntegralAbstract implements UnidirectionalEditDistance {

  private static final boolean DEBUG = false;

  protected int[] mWorkspace = new int[ActionsHelper.ACTIONS_START_INDEX];

  static final int LEFT_ARM = 0;
  static final int RIGHT_ARM = 1;

  static final int OVERLAP_GAP = RealignParamsImplementation.CG_OVERLAP;
  static final int SMALL_GAP = RealignParamsImplementation.CG_SMALL_GAP;
  static final int LARGE_GAP = RealignParamsImplementation.CG_LARGE_GAP;
  static final int OVERLAP2_GAP = RealignParamsImplementation.CG_OVERLAP2;

  /** This is indexed by the above three constants. */
  static final String[] GAP_NAME = {"Overlap", "Small", "Large", "Overlap2"};

  /** width of each number in the <code>toString</code> output */
  static final int FIELD_WIDTH = 9;

  protected final RealignParams mParams;

  protected final double mMatchProb;
  protected final double mMisMatchProb;
  protected final double mDeleteOpenProb;
  protected final double mInsertOpenProb;
  protected final double mDeleteExtendProb;
  protected final double mInsertExtendProb;
  protected final double mOneMinusDeleteExtend;
  protected final double mOneMinusDeleteOpen;
  protected final double mOneMinusInsertExtend;
  protected final double mOneMinusDeleteInsertOpen;
  protected final int[] mGapStart;
  protected final double[][] mGapProb;

  protected Environment mEnv;

  protected double[][] mInsert;
  protected double[][] mMatch;
  protected double[][] mDelete;

  /** The height of the matrix less 1 (so valid rows are from 0 .. <code>mLength</code> inclusive) */
  protected int mLength;

  /** The minimum width of the matrix (i.e., the width of the first row) */
  protected final int mWidth;

  /** Maximum variation in gap size, for recording statistics. */
  private static final int MAX_GAP = 7;
  /** Which arm of statistics this read should count towards. */
  private int mStatsArm;
  private final int[][] mStats;

  private byte[] mRead;

  /**
   * Index into <code>mRowOffsets</code> to determine which arm we are processing.
   * For 35 bp reads this means:
   * 0 means left arm <code>(5-overlap-10-smallgap-10-largegap-10)</code>.
   * 1 means right arm <code>(10-largegap-10-smallgap-10-overlap-5)</code>.
   * For 29 bp reads this means:
   * 0 means left arm <code>(10-overlap2-19)</code>.
   * 1 means right arm <code>(10-overlap2-19)</code>.
   */
  private int mArm;

  /**
   * Template offset of start of each row (relative to expected start of read)
   * Indexed by the arm, then by the one-based read position.
   */
  private final int[][] mRowOffsets;

  /**
   * <code>mGap[arm][pos]</code> gives the CG gap that happens just before <code>pos</code>
   * (where <code>pos</code> is a one-based read position).
   * The resulting value will be -1 (no gap), <code>OVERLAP_GAP</code>, <code>SMALL_GAP</code>, <code>LARGE_GAP</code>, or <code>OVERLAP2_GAP</code>.
   */
  private final int[][] mGap;

  private final int mUnknownsPenalty;

  /** The default size of the insert region for CG data */
  public static final int CG_INSERT_REGION_DEFAULT_SIZE = 5;

  /**
   * @param maxShift maximum distance that start and end of alignment can move.
   * @param params CG probabilities.
   * @param unknownsPenalty penalty for an unknown nucleotide - can only either be 0 or 1 for CG.
   */
  public CgGotohEditDistance(final int maxShift, final RealignParams params, int unknownsPenalty) {
    this(maxShift, params, unknownsPenalty, false);
  }

  /**
   * @param maxShift maximum distance that start and end of alignment can move.
   * @param params CG probabilities.
   * @param unknownsPenalty penalty for an unknown nucleotide - can only either be 0 or 1 for CG.
   * @param v2 true for CG version 2 read structure, otherwise assume version 1 read structure
   */
  public CgGotohEditDistance(final int maxShift, final RealignParams params, int unknownsPenalty, boolean v2) {
    assert params.completeGenomics();
    mEnv = null; // set by each edit distance call.
    mLength = v2 ? CgUtils.CG2_RAW_READ_LENGTH : CgUtils.CG_RAW_READ_LENGTH;
    mWidth = maxShift * 2 + 1;
    mParams = params;
    mMatchProb = Math.exp(mParams.matchLn());
    mMisMatchProb = Math.exp(mParams.misMatchLn());
    mDeleteOpenProb = Math.exp(mParams.deleteOpenLn());
    mInsertOpenProb = Math.exp(mParams.insertOpenLn());
    mDeleteExtendProb = Math.exp(mParams.deleteExtendLn());
    mInsertExtendProb = Math.exp(mParams.insertExtendLn());
    mOneMinusDeleteExtend = 1.0 - mDeleteExtendProb;
    mOneMinusDeleteOpen = 1.0 - mDeleteOpenProb;
    mOneMinusInsertExtend = 1.0 - mInsertExtendProb;
    mOneMinusDeleteInsertOpen = 1.0 - mDeleteOpenProb - mInsertOpenProb;
    mGapStart = new int[GAP_NAME.length];
    mGapProb = new double[GAP_NAME.length][];
    for (int gap = 0; gap < GAP_NAME.length; gap++) {
      mGapStart[gap] = mParams.gapStart(gap);
      mGapProb[gap] = new double[mParams.gapEnd(gap) - mParams.gapStart(gap) + 1];
      for (int i = 0; i < mGapProb[gap].length; i++) {
        mGapProb[gap][i] = Math.exp(mParams.gapFreqLn(gap, mGapStart[gap] + i));
      }
    }
    mMatch = makeDouble(mLength + 1);
    mInsert = makeDouble(mLength + 1);
    mDelete = makeDouble(mLength + 1);

    // set up CG gap positions
    mGap = new int[2][];
    mGap[LEFT_ARM] = new int[mLength + 1];
    Arrays.fill(mGap[LEFT_ARM], -1);
    mGap[RIGHT_ARM] = new int[mLength + 1];
    Arrays.fill(mGap[RIGHT_ARM], -1);
    if (!v2) {
      mGap[LEFT_ARM][6] = OVERLAP_GAP;
      mGap[LEFT_ARM][16] = SMALL_GAP;
      mGap[LEFT_ARM][26] = LARGE_GAP;
      mGap[RIGHT_ARM][11] = LARGE_GAP;
      mGap[RIGHT_ARM][21] = SMALL_GAP;
      mGap[RIGHT_ARM][31] = OVERLAP_GAP;
    } else {
      mGap[LEFT_ARM][11] = OVERLAP2_GAP;
      mGap[RIGHT_ARM][11] = OVERLAP2_GAP;
    }
    // set up CG row offsets for the most common gap sizes
    mRowOffsets = new int[2][];
    mRowOffsets[LEFT_ARM] = makeRowOffsets(LEFT_ARM);
    mRowOffsets[RIGHT_ARM] = makeRowOffsets(RIGHT_ARM);

    mStats = new int[2 * GAP_NAME.length][];
    for (int i = 0; i < mStats.length; i++) {
      mStats[i] = new int[MAX_GAP];
    }
    mUnknownsPenalty = unknownsPenalty == 0 ? 0 : 1; //cg ed only supports 1/1/1/1 or 1/1/1/0 penalties, really.
  }

  /**
   * Change the environment and recalculate all alignment paths.
   *
   * @param env the new read and template information.
   * @param leftArm true for left arm, false for right arm.
   */
  public void setEnv(final Environment env, final boolean leftArm) {
    assert env.readLength() == mLength;
    mEnv = env;
    mArm = leftArm ? LEFT_ARM : RIGHT_ARM;
    calculateProbabilities();
    assert globalIntegrity();
  }

  /**
   * Calculate all the probabilities in the matrix.
   */
  protected void calculateProbabilities() {
    assert mLength == CgUtils.CG_RAW_READ_LENGTH || mLength == CgUtils.CG2_RAW_READ_LENGTH : "mLength=" + mLength; // it should be a CG read
    calculateInitialRow(0, mDeleteOpenProb / mWidth, (1.0 - mDeleteOpenProb) / mWidth);
    for (int row = 1; row <= mLength; row++) {
      final int gap = mGap[mArm][row];
      if (gap < 0) {
        calculateRow(row); // normal boring row
      } else {
        calculateCGGap(row, gap); // an exciting CG gap between row-1 and row
      }
    }
  }

  /**
   * How far a given row is shifted along the template, relative to
   * the expected start position of the read.  For CG, this accounts for
   * the CG gaps - it is used to keep the middle of the alignment band
   * in line with the most common offsets of the CG gaps.
   *
   * @param row one-based read position.
   * @return an offset along the template, from where the read is expected to start.
   *         Smallest value will be for row 0.
   *         Largest value will be for the last row.
   *         But the values for intermediate rows may move backwards
   *         and forwards between those bounds.
   */
  protected int rowOffset(final int row) {
    return mRowOffsets[mArm][row];
  }

  int[] makeRowOffsets(final int arm) {
    final int maxShift = mWidth >> 1;
      final int[] result = new int[mLength + 1];
      int offset = -(maxShift + 1);
      for (int row = 0; row <= mLength; row++) {
        final int gap = mGap[arm][row];
        if (gap == OVERLAP_GAP) {
          offset -= 2;
        } else if (gap == OVERLAP2_GAP) {
          offset -= 3;
        } else if (gap == LARGE_GAP) {
          offset += 6;
        }
        result[row] = offset;
        offset++;
      }
      return result;
  }

  protected double calculateDelete(final int i, final int j) {
    if (j == mWidth) {
      return 0.0;
    }
    final double sum = mDeleteExtendProb * mDelete[i][j] + mDeleteOpenProb * mMatch[i][j];
    return sum / 4.0;
  }

  /**
   * This is similar to <code>calculateDelete</code>, but always uses the
   * open-delete penalty, never the extend-delete penalty.
   * It should be used whenever a CG gap is crossed, even if the gap is zero,
   * because the chemistry is discontinuous at each gap.
   * @param i row
   * @param j column
   * @return probability.
   */
  protected double calculateOpenDelete(final int i, final int j) {
    return (mDelete[i][j] + mMatch[i][j]) * mDeleteOpenProb / 4.0;
  }

  /**
   * Calculate the forward diagonal probabilities, starting from cell <code>(i,j)</code>.
   * This does NOT include the match/mismatch probability of the destination cell.
   *
   * @param i zero-based read position
   * @param j column
   * @return probability
   */
  protected double calculateMatch(final int i, final int j) {
    final double del = mDelete[i][j] * mOneMinusDeleteExtend;
    final double match = mMatch[i][j] * mOneMinusDeleteInsertOpen;
    final double ins = mInsert[i][j] * mOneMinusInsertExtend;
    return ins + del + match;
  }

  // similar to calculateMatch, but always uses the delete open penalty.
  protected double calculateMatchCGGap(final int i, final int j) {
    final double del = mDelete[i][j] * mOneMinusDeleteOpen;
    final double match = mMatch[i][j] * mOneMinusDeleteInsertOpen;
    final double ins = mInsert[i][j] * mOneMinusInsertExtend;
    return del + match + ins;
  }

  protected double calculateInsert(final int i, final int j) {
    if (i == mLength || j < 0) {
      return 0.0;
    }
    return mInsertExtendProb * mInsert[i][j] + mInsertOpenProb * mMatch[i][j];
  }

  protected void calculateRow(final int i) {
    for (int j = 0; j < mWidth; j++) {
      mDelete[i][j] = calculateDelete(i - 1, j + 1);
      mMatch[i][j] = calculateMatch(i - 1, j) * matchEq(i - 1, j);
      mInsert[i][j] = calculateInsert(i, j - 1);
    }
  }

  /**
   * Calculates a CG row after a variable size gap.
   * @param row the first row after the gap.
   * @param whichGap which gap to handle (see <code>RealignParams</code>).
   */
  protected void calculateCGGap(final int row, final int whichGap) {
    final int gapStart = mGapStart[whichGap];
    final double[] gapProb = mGapProb[whichGap];
    final int gapEnd = gapStart + gapProb.length; // one PAST the final gap
    //    final byte re = mEnv.read(row - 1);
    final int lastCol = mWidth - 1;
    final int thisRowOffset = rowOffset(row);
    final int offset = thisRowOffset - rowOffset(row - 1);

    for (int j = 0; j <= lastCol; j++) {
      // delete: we sum the probabilities over all gap sizes
      double del = 0.0;
      for (int gapSize = gapStart; gapSize < gapEnd; gapSize++) {
        final int prevCol = j - (gapSize - offset);
        if (0 <= prevCol && prevCol < mWidth) {
          final double from = calculateOpenDelete(row - 1, prevCol);
          del += from * gapProb[gapSize - gapStart];
        }
      }
      mDelete[row][j] = del;

      // match or mismatch: we sum the probabilities over all gap sizes
      double mm = 0.0;
      for (int gapSize = gapStart; gapSize < gapEnd; gapSize++) {
        final int prevCol = j - (gapSize - offset) - 1;
        //        final byte te = mEnv.template(thisRowOffset + j);
        if (0 <= prevCol && prevCol < mWidth) {
          mm += calculateMatchCGGap(row - 1, prevCol) * gapProb[gapSize - gapStart];
        }
      }
      mMatch[row][j] = mm * matchEq(row - 1, j);

      // insert is the same as usual
      mInsert[row][j] = calculateInsert(row, j - 1);
    }
  }

  protected void calculateInitialRow(final int initRow, final double delete, final double match) {
    for (int j = 0; j < mWidth; j++) {
      mDelete[initRow][j] = delete;
      mMatch[initRow][j] = match;
      mInsert[initRow][j] = 0.0;
    }
  }

  /**
   * Calculate the match/mismatch cost at a given cell of the matrix.
   * @param i zero-based read position
   * @param j zero-based column position (<code>0 .. mWidth - 1</code>)
   * @return natural log of the probability
   */
  protected double matchEq(final int i, final int j) {
    final int templatePos = rowOffset(i + 1) + j;
    final byte te = mEnv.template(templatePos);
    final int absTemplatePos = mEnv.absoluteTemplatePosition(templatePos);
    return matchEqTe(i, te, absTemplatePos >= 0 && absTemplatePos < mEnv.templateLength());
  }

  /**
   * Calculate the match/mismatch cost at a given row and for a specified nucleotide.
   * @param i zero-based read position
   * @param nt nucleotide (assumed to be on template, 0=N...4=T).
   * @param isWithinTemplate true if the position is on the template
   * @return the probability
   */
  double matchEqTe(int i, byte nt, boolean isWithinTemplate) {
    if (!isWithinTemplate) {
      return 0.25;
    }
    final byte re = mEnv.read(i);
    final double q = mEnv.quality(i);
    final double q3 = q / 3.0;
    //System.err.println("         i=" + i + " j=" + j + " te=" + te + " re=" + re + " sum=" + sum);
    final double incr;
    if (nt == 0 || re == 0) {
      if (mUnknownsPenalty > 0) {
        incr = 0.25;
      } else {
        incr = mMatchProb * (1.0 - q) + mMisMatchProb * q3;
      }
    } else if (nt == re) {
      incr = mMatchProb * (1.0 - q) + mMisMatchProb * q3;
      //System.err.println("match    i=" + i + " j=" + j + " incr=" + incr + " res=" + (sum + incr));
    } else {
      incr = mMatchProb * q3 + mMisMatchProb * (1.0 - q3) / 3.0;
      //System.err.println("mismatch i=" + i + " j=" + j + " x=" + x + " y=" + y + " incr=" + incr + " res=" + (sum + incr));
    }
    return incr;
  }

  /**
   * Sharpen creating a two dimensional array.
   * @param rows length of first index.
   * @return the new array.
   */
  protected final double[][] makeDouble(final int rows) {
    final double[][] res = new double[rows][];
    for (int i = 0; i < rows; i++) {
      res[i] = new double[mWidth];
    }
    return res;
  }

  /**
   * Pads out two string by inserting spaces between them.
   *
   * @param first left string
   * @param length total length of the desired output string
   * @param last right string
   * @return the padded string
   */
  public static String extend(final String first, final int length, final String last) {
    final int len = first.length() + last.length();
    if (len < length) {
      return first + StringUtils.getSpaceString(length - len) + last;
    } else {
      return first + last;
    }
  }

  // formats like 1.234e-04
  static String format(final double x) {
    final String result;
    final int expWidth = 4;
    final int dp = (FIELD_WIDTH - expWidth) / 2;
    if (x == 0.0) {
      // avoid outputting -0.0
      result = "";
    } else if (Double.isInfinite(x)) {
      result = (x < 0.0 ? "-" : "") + "inf";
    } else {
      final String fmt = "%1$01." + dp + "e";
      result = String.format(Locale.ROOT, fmt, x);
    }
    return extend("", FIELD_WIDTH, result);
  }

  @Override
  public final void toString(final StringBuilder sb) {
    final char[] dnaChars = DNA.valueChars();
    final int rowStart = rowOffset(0);
    final int rowEnd = rowOffset(mLength) + mWidth;
    printTemplateRow(sb, rowStart, rowEnd, mEnv, "ScoreMatrix");
    for (int row = 0; row <= mLength; row++) {
      if (mGap[mArm][row] >= 0) {
        // show a horizontal line where the CG gap/overlap happens.
        sb.append("      ");
        for (int i = 0; i < (FIELD_WIDTH * 3 + 1) * (rowOffset(row) - rowStart + mWidth); i++) {
          sb.append('-');
        }
        sb.append(LS);
      }
      sb.append("[").append(extend("", 3, row + "")).append("]");
      final char re = row == 0 ? ' ' : dnaChars[mEnv.read(row - 1)];
      sb.append(re);
      for (int tPos = 0; tPos < rowOffset(row) - rowStart + mWidth; tPos++) {
        final String sep = row > 0 && (tPos + rowStart + 1) % 5 == 0 ? ("" + re) : "|";
        if (tPos < rowOffset(row) - rowStart) {
          // indent the row, so we line up with the template
          sb.append(extend("", 3 * FIELD_WIDTH + 1, sep));
        } else {
          final int col = tPos - (rowOffset(row) - rowStart);
          sb.append(format(mInsert[row][col]));
          sb.append(format(mMatch[row][col]));
          sb.append(format(mDelete[row][col]));
          sb.append(sep);
        }
      }
      sb.append(LS);
    }
    printTemplateRow(sb, rowStart, rowEnd, mEnv, "");
  }

  static void printTemplateRow(final StringBuilder sb, final int rowStart, final int rowEnd, final Environment env, final String msg) {
    // print the header row, showing the template.
    final char[] dnaChars = DNA.valueChars();
    sb.append(extend(msg, 7 + 3 * FIELD_WIDTH, "|"));
    for (int pos = rowStart + 1; pos <= rowEnd; pos++) {
      sb.append(dnaChars[env.template(pos)]);
      sb.append(extend("", 2 * FIELD_WIDTH - 1, Integer.toString(env.absoluteTemplatePosition(pos))));
      sb.append(extend("", FIELD_WIDTH + 1, "|"));
    }
    sb.append(LS);
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(0 < mLength);
    Exam.assertTrue(0 < mWidth);
    Exam.assertTrue(rowOffset(0) < rowOffset(mLength));
    Exam.assertEquals(mLength + 1, mMatch.length);
    Exam.assertEquals(mLength + 1, mDelete.length);
    Exam.assertEquals(mLength + 1, mInsert.length);
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (int i = 0; i < mLength; i++) {
      // the three matrices have exactly the same shape.
      Exam.assertEquals(mWidth, mMatch[i].length);
      Exam.assertEquals(mWidth, mDelete[i].length);
      Exam.assertEquals(mWidth, mInsert[i].length);
      for (int j = 0; j < mWidth; j++) {
        final double vi = mInsert[i][j];
        Exam.assertTrue("i=" + i + " j=" + j + " vi=" + vi, 0.0 <= vi && vi < 1.0 && !Double.isNaN(vi));
        final double vm = mMatch[i][j];
        Exam.assertTrue("i=" + i + " j=" + j + " vm=" + vm, 0.0 <= vm && vm <= 1.0 && !Double.isNaN(vm));
        final double vd = mDelete[i][j];
        Exam.assertTrue("i=" + i + " j=" + j + " vd=" + vd, 0.0 <= vd && vd <= 1.0 && !Double.isNaN(vd));
      }
    }
    return true;
  }

  /**
   * @param readEnd the last position in the read (0 based exclusive)
   * @param reverse if true then reverse the actions array at the end
   * @param workSpace0 where to store the actions array.
   * @return the actions array - in some cases will be <code>workSpace0</code> in others it will be newly created.
   */
  public int[] goBackwards(final int readEnd, final boolean reverse, final int[] workSpace0) {
    int bestPos = 0;
    double bestScore = 0.0;
    for (int i = 0; i < mWidth; i++) {
      final double tempScore = mMatch[mLength][i] + mDelete[mLength][i];
      if (tempScore > bestScore) {
        bestScore = tempScore;
        bestPos = i;
      }
    }
    final int endAlignment = bestPos;
    final int workLen = ActionsHelper.ACTIONS_START_INDEX + 1 + ((readEnd + mWidth) * 2) / ActionsHelper.ACTIONS_PER_INT;
    final int[] workSpace;
    if (workSpace0.length < workLen) {
      workSpace = new int[workLen];
    } else {
      workSpace = workSpace0;
    }
    ActionsHelper.clear(workSpace);
    int readPos = readEnd;
    int refPos = endAlignment;
    int previousAction = -1;
    while (readPos > 0) {
      final double dd = mDelete[readPos][refPos];
      final double mm = mMatch[readPos][refPos];
      final double ii = mInsert[readPos][refPos];
      final double best = Math.max(dd, Math.max(mm, ii));
      final int thisOffset = rowOffset(readPos);
      final int nextOffset = rowOffset(readPos - 1); // next row we are going to
      //System.err.println("(" + readPos + "," + refPos + ") rowOffset=" + thisOffset + " nextOffset=" + nextOffset + " as=" + mWorkspace[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
      if (ii == best) {
        // NOTE: ii is zero when refPos==0, so we should never make refPos negative here.
        ActionsHelper.prepend(workSpace, 1, ActionsHelper.DELETION_FROM_REFERENCE, previousAction == ActionsHelper.DELETION_FROM_REFERENCE ? 1 : 2);
        previousAction = ActionsHelper.DELETION_FROM_REFERENCE;
        refPos--;
        assert refPos >= 0;
      } else if (dd == best) {
        ActionsHelper.prepend(workSpace, 1, ActionsHelper.INSERTION_INTO_REFERENCE, previousAction == ActionsHelper.INSERTION_INTO_REFERENCE ? 1 : 2);
        previousAction = ActionsHelper.INSERTION_INTO_REFERENCE;
        if (mGap[mArm][readPos] >= 0) {
          final int bestGap = findBestGap(readPos, refPos, 0);
          // System.err.println("readPos=" + readPos + " delete bestGap=" + bestGap);
          if (bestGap != 0) {
            final int cgCmd = bestGap < 0 ? ActionsHelper.CG_OVERLAP_IN_READ : ActionsHelper.CG_GAP_IN_READ;
            ActionsHelper.prepend(workSpace, Math.abs(bestGap), cgCmd, 0);
          }
          refPos = refPos + thisOffset - bestGap - nextOffset;
          assert 0 <= refPos;
          assert refPos < mWidth;
        } else {
          // NOTE: dd is zero when refPos==mWidth-1, so refPos should stay < mWidth.
          refPos++;
          assert refPos < mWidth;
        }
        readPos--;
      } else {
        final int templatePos = thisOffset + refPos;
        final int cmd;
        int penalty = 0;
        if (mEnv.absoluteTemplatePosition(templatePos) < 0 || mEnv.absoluteTemplatePosition(templatePos) >= mEnv.templateLength()) {
          cmd = ActionsHelper.MISMATCH;
        } else {
          final byte te = mEnv.template(templatePos);
          final byte re = mEnv.read(readPos - 1);
          if (re == 0) {  //check read for N first, so it doesn't need to be evaluated for read delta
            cmd = ActionsHelper.UNKNOWN_READ;
            penalty = mUnknownsPenalty;
          } else if (te == 0) {
            cmd = ActionsHelper.UNKNOWN_TEMPLATE;
            penalty = mUnknownsPenalty;
          } else {
            cmd = (te == re) ? ActionsHelper.SAME : ActionsHelper.MISMATCH;
          }
        }
        ActionsHelper.prepend(workSpace, 1, cmd, cmd == ActionsHelper.MISMATCH ? 1 : penalty);
        if (mGap[mArm][readPos] >= 0) {
          final int bestGap = findBestGap(readPos, refPos, 1);
          //System.err.println("readPos=" + readPos + " equality bestGap=" + bestGap);
          if (bestGap != 0) {
            final int cgCmd = bestGap < 0 ? ActionsHelper.CG_OVERLAP_IN_READ : ActionsHelper.CG_GAP_IN_READ;
            ActionsHelper.prepend(workSpace, Math.abs(bestGap), cgCmd, 0);
          }
          refPos = refPos + thisOffset - 1 - bestGap - nextOffset;
          assert 0 <= refPos;
          assert refPos < mWidth;
        }
        readPos--;
        previousAction = cmd;
      }
    }
    if (DEBUG) {
      System.err.println("start := " +  mEnv.absoluteTemplatePosition(rowOffset(0) + refPos + 1));
      System.err.println("read=" + Arrays.toString(mRead));
      System.err.print("tmpl=[");
      for (int i = 0; i < 40; i++) {
        System.err.print(mEnv.template(i) + ", ");
      }
      System.err.println("]");
    }
    ActionsHelper.setZeroBasedTemplateStart(workSpace, mEnv.absoluteTemplatePosition(rowOffset(0) + refPos + 1));
    if (reverse) {
      // reverse the alignment actions and fix the start position.
      ActionsHelper.reverse(workSpace);
      ActionsHelper.setZeroBasedTemplateStart(workSpace, mEnv.absoluteTemplatePosition(rowOffset(mLength) + bestPos));
    }
    return workSpace;
  }

  /**
   * @param row one-based read position
   * @param col equivalent to template position minus row offset.
   * @param match 1 to look for the best match gap, 0 to look for the best delete gap.
   * @return the size of the most probable gap (negative means overlap).
   */
  private int findBestGap(int row, int col, int match) {
    final int whichGap = mGap[mArm][row];
    assert whichGap >= 0;
    final int gapStart = mGapStart[whichGap];
    final double[] gapProb = mGapProb[whichGap];
    final int gapEnd = gapStart + gapProb.length; // one PAST the final gap
    final int thisRowOffset = rowOffset(row);
    final int offset = thisRowOffset - rowOffset(row - 1);
    int bestGap = 0;
    double bestProb = 0.0;
    for (int gapSize = gapStart; gapSize < gapEnd; gapSize++) {
      final int prevCol = col - (gapSize - offset) - match;
      if (0 <= prevCol && prevCol < mWidth) {
        final double prob = ((match == 1)
            ? calculateMatchCGGap(row - 1, prevCol)
                : calculateOpenDelete(row - 1, prevCol))
                * gapProb[gapSize - gapStart];
        if (prob > bestProb) {
          bestGap = gapSize;
          bestProb = prob;
        }
      }
    }
    assert bestProb > 0.0;
    mStats[mStatsArm * GAP_NAME.length + whichGap][bestGap - gapStart]++;
    return bestGap;
  }

  @Override
  public void logStats() {
    final FormatReal fmt = new FormatReal(1, 6);
    final StringBuilder sb = new StringBuilder();
    sb.append("CgGotohEditDistance mismatch=");
    sb.append(fmt.format(mMisMatchProb));
    sb.append(" delopen=");
    sb.append(fmt.format(mDeleteOpenProb));
    sb.append(" delextend=");
    sb.append(fmt.format(mDeleteExtendProb));
    sb.append(" insopen=");
    sb.append(fmt.format(mInsertOpenProb));
    sb.append(" insextend=");
    sb.append(fmt.format(mInsertExtendProb));
    sb.append(StringUtils.LS);
    sb.append("Gap Histograms").append(StringUtils.LS);
    for (int arm = 0; arm < 2; arm++) {
      final int startGap = 0;
      final int endGap = 3;
      for (int whichGap = startGap; whichGap < endGap; whichGap++) {
        final int index = arm * GAP_NAME.length + whichGap;
        sb.append(arm == LEFT_ARM ? "Left" : "Right").append(GAP_NAME[whichGap]).append(":\t");
        final int gapStart = mGapStart[whichGap];
        final int gapEnd = gapStart + mGapProb[whichGap].length;
        sb.append(gapStart).append("..").append(gapEnd - 1);
        for (int width = 0; width < gapEnd - gapStart; width++) {
          sb.append("\t").append(mStats[index][width]);
        }
        sb.append(StringUtils.LS);
      }
    }
    Diagnostic.developerLog(sb.toString());
  }

  @Override
  public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
    // TODO: get read qualities from somewhere
    assert read.length == mLength;
//    System.err.println("aligning " + rlen + "@" + zeroBasedStart + " : " + cgLeft);
    final EnvironmentImplementation env;
    final boolean isInverted = !cgLeft && rlen == CgUtils.CG_RAW_READ_LENGTH;
    mStatsArm = cgLeft ? LEFT_ARM : RIGHT_ARM;
    if (isInverted) {
      // we reverse the read and the template, and then pretend that this is a left arm read.
      // we add an offset because the most common CG alignment size along the template is 35 - 2 + 6 = 39.
      env = new InvertedEnvironment(new EnvironmentImplementation(maxShift, template, zeroBasedStart + CgUtils.CG_EXPECTED_LENGTH_OFFSET, read, null));
    } else {
      env = new EnvironmentImplementation(maxShift, template, zeroBasedStart, read, null);
    }
    setEnv(env, true);
    mRead = read;
    mWorkspace = goBackwards(mLength, isInverted, mWorkspace);
//    System.err.println("workspace says aligns @" + mWorkspace[ActionsHelper.TEMPLATE_START_INDEX] + " score=" + mWorkspace[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    return mWorkspace;
  }

  @Override
  public int[] calculateEditDistanceFixedBoth(byte[] read, int readStartPos, int readEndPos,
      byte[] template, int templateStart, int templateEnd, int maxScore, int maxShift) {
    return null;
  }

  @Override
  public int[] calculateEditDistanceFixedEnd(byte[] read, int readStartPos, int readEndPos, byte[] template,
      int templateExpectedStartPos, int templateEndPos, int maxScore, int maxShift) {
    return null;
  }

  @Override
  public int[] calculateEditDistanceFixedStart(byte[] read, int readStartPos, int readEndPos,
      byte[] template, int templateStartPos, int maxScore, int maxShift) {
    return null;
  }
}
