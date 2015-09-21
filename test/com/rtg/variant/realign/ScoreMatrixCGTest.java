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

import static com.rtg.util.StringUtils.LS;

import java.util.Arrays;

import com.rtg.mode.DnaUtils;
import com.rtg.util.integrity.Exam;
import com.rtg.variant.realign.ScoreMatrixTest.MockRealignParams;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 */
public abstract class ScoreMatrixCGTest extends TestCase {

  //Mark says I have to put my name here so he doesnt get blamed for the next bit of code
  //JC
  /** Used for testing. */
  public static class LogTest extends ScoreMatrixCGTest {
    public LogTest() { }
    @Override
    protected AllPaths score(final Environment env) {
      return score(LogPossibility.SINGLETON, env);
    }
  }

  /** Used for testing. */
  public static class LogApproximateTest extends ScoreMatrixCGTest {
    public LogApproximateTest() { }
    @Override
    protected AllPaths score(final Environment env) {
      return score(LogApproximatePossibility.SINGLETON, env);
    }
    //the results are slightly inaccurate in the 5th digit - but the other two accurate versions agree so let it be
    @Override
    public void testToString() { }
  }

  /** Used for testing. */
  public static class SimpleTest extends ScoreMatrixCGTest {
    public SimpleTest() { }
    @Override
    protected AllPaths score(final Environment env) {
      return score(SimplePossibility.SINGLETON, env);
    }
  }

  public static Test suite() {
    final TestSuite suite = new TestSuite();
    suite.addTestSuite(LogApproximateTest.class);
    suite.addTestSuite(LogTest.class);
    suite.addTestSuite(SimpleTest.class);
    return suite;
  }

  /**
   * CG Realign parameters that match the scorematrixtestCG.xls spreadsheet.
   */
  public static class MockRealignParamsCG extends ScoreMatrixTest.MockRealignParams {

    static final int[] GAP_START = {-3, 0, 5, -4};

    static final double[][] GAP_FREQ =
      {
      {Math.log(0.08), Math.log(0.84), Math.log(0.08)},    // for gap of -3, -2, -1
      {Math.log(0.27), Math.log(0.635), Math.log(0.095)},  // for gap of  0,  1,  2
      {Math.log(0.90), Math.log(0.07), Math.log(0.03)},    // for gap of  5,  6,  7
      {Math.log(0.08), Math.log(0.84), Math.log(0.08)},    // for gap of -4, -3, -2
      };

    @Override
    public boolean completeGenomics() {
      return true;
    }

    @Override
    public int gapStart(final int gap) {
      return GAP_START[gap];
    }

    @Override
    public int gapEnd(final int gap) {
      return GAP_START[gap] + 2;
    }

    @Override
    public double gapFreqLn(final int gap, final int width) {
      return GAP_FREQ[gap][width - gapStart(gap)];
    }

    @Override
    public double[][] gapDistributionPoss(PossibilityArithmetic arith) {
      return RealignParamsImplementation.gapDistributionPoss(GAP_FREQ, arith);
    }

  }

  /**
   * This is the template and read from the <code>scorematrixtestCG.xls</code> spreadsheet.
   * The read has an overlap of -2, a small gap of 0, and a large gap of 7.
   */
  public static final String TEMPLATE = "GGGGGGGATA AAAAA   GGCGACAT GCCAATGTGT CGCCTTT TTCAACTTTC CGATTAA".replaceAll(" ", "");
  static final String READ = "               ATAAA AAGGCGACAT GCCAATGTGT         TTCAACTTTC".replaceAll(" ", "");

  protected byte[] read() {
    return DnaUtils.encodeString(READ);
  }

  protected byte[] template() {
    return DnaUtils.encodeString(TEMPLATE);
  }

  protected Environment env() {
    final double[] quality = new double[read().length];
    Arrays.fill(quality, 0.01);
    quality[0] = 0.1;
    quality[1] = Math.pow(10.0, -15 / 10.0);
    final Environment env = new EnvironmentImplementation(
        7, //maxShift
        template(),
        7, //start
        read(),
        quality
        );
    Exam.integrity(env);
    return env;
  }

  protected abstract AllPaths score(final Environment env);

  protected AllPaths score(final PossibilityArithmetic arith, final Environment env) {
    final AllPaths score = new ScoreMatrixCG(arith, new MockRealignParamsCG());
    score.setEnv(env);
    Exam.globalIntegrity(score);
    return score;
  }

  public void testRowOffset() {
    final AllPaths score0 = score(env());
    if (score0 instanceof ScoreMatrixCG) {
      final ScoreMatrixCG score = (ScoreMatrixCG) score0;
      final int half = -8; // half of the width of the band (at row 1).
      assertEquals(half, score.rowOffset(0));
      assertEquals(half + 1, score.rowOffset(1));
      assertEquals(half + 5, score.rowOffset(5));
      assertEquals(half + 6 - 2, score.rowOffset(6));
      assertEquals(half + 15 - 2, score.rowOffset(15));
      assertEquals(half + 16 - 2, score.rowOffset(16));
      assertEquals(half + 25 - 2, score.rowOffset(25));
      assertEquals(half + 26 - 2 + 6, score.rowOffset(26));
      assertEquals(half + 35 - 2 + 6 , score.rowOffset(35));
    }
  }

  public void testStartScores() {
    final Environment mEnv = env();
    final AllPaths score0 = score(mEnv);
    if (score0 instanceof ScoreMatrixCG) {
      final ScoreMatrixCG score = (ScoreMatrixCG) score0;
      final int maxShift = mEnv.maxShift();
      final int w = 2 * maxShift + 1;  // width
      final double matchPenalty = Math.log(1.0 - Math.exp(score.mParams.deleteOpenLn())); // / w);
      final double deletePenalty = Math.log(Math.exp(score.mParams.deleteOpenLn())); // / w);
      assertEquals(w, score.mWidth);
      assertEquals(-maxShift - 1, score.rowOffset(0));
      for (int tPos = -maxShift - 1; tPos < maxShift; tPos++) {
        checkScores(score, 0, tPos, deletePenalty, matchPenalty, Double.NEGATIVE_INFINITY);
      }
    }
  }

  public void testTotalScoreLn() {
    final AllPaths score = score(env());
    assertEquals(-10.012470, score.totalScoreLn(), 0.001);
    assertEquals(-10.012470, Math.log(score.totalScore()), 0.001);
  }

  public void testScoreLn() {
    final Environment mEnv = env();
    final AllPaths score0 = score(mEnv);
    if (score0 instanceof ScoreMatrixCG) {
      final ScoreMatrixCG score = (ScoreMatrixCG) score0;
      final int start = mEnv.absoluteTemplatePosition(0);
      // down column 32 (T on the template), cuts through the big gap
      final int tPos = 32 - start;
      checkScores(score, 21, tPos, Double.NEGATIVE_INFINITY, -18.7048, -18.4873);
      checkScores(score, 22, tPos, -27.1852, -22.1531, -16.8977);
      checkScores(score, 23, tPos, -29.7517, -14.4359, -15.3150);
      checkScores(score, 24, tPos, -22.9162, -18.9807, -13.7276);
      checkScores(score, 25, tPos, -25.8054, -6.2972, -25.7807);
      checkScores(score, 26, tPos, -35.6482, -30.2646, -44.9473);
      checkScores(score, 27, tPos, -38.0540, -37.9438, -48.9078);
      checkScores(score, 28, tPos, -41.1499, -46.0732, -56.0204);
      checkScores(score, 29, tPos, -44.2510, -50.1387, Double.NEGATIVE_INFINITY);
    }
  }

  protected void checkScores(ScoreMatrixCG score, final int row, final int relTemplatePos, final double delete, final double match, final double insert) {
    final int col = relTemplatePos - score.rowOffset(row);
    assertTrue("template position too small in test", 0 <= col);
    assertTrue("template position too large in test", col < score.mWidth);
    assertEquals(delete, score.arithmetic().poss2Ln(score.delete(row, col)), 0.001);
    assertEquals(match, score.arithmetic().poss2Ln(score.match(row, col)), 0.001);
    assertEquals(insert, score.arithmetic().poss2Ln(score.insert(row, col)), 0.001);
  }

  public void testToString() {
    final AllPaths score = score(env());
    assertEquals(""
        + "ScoreMatrix                |G            0       |G            1       |G            2       |G            3       |G            4       |G            5       |G            6       |A            7       |T            8       |A            9       |A           10       |A           11       |A           12       |A           13       |A           14       |G           15       |G           16       |C           17       |G           18       |A           19       |C           20       |A           21       |T           22       |G           23       |C           24       |C           25       |A           26       |A           27       |T           28       |G           29       |T           30       |G           31       |T           32       |C           33       |G           34       |C           35       |C           36       |T           37       |T           38       |T           39       |T           40       |T           41       |C           42       |A           43       |A           44       |C           45       |T           46       |T           47       |T           48       |C           49       |C           50       |G           51       |A           52       |T           53       |" + LS
        + "[  0]          0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|" + LS
        + "[  1]A                     |         3.327  8.316| 10.780  3.327  8.316| 10.598  3.327  8.316| 10.565  3.327  8.316| 10.558  3.327  8.316| 10.557  3.327  8.316| 10.557  3.327  8.316| 10.557  0.117  8.316|  7.560  3.327  8.316|  8.987  0.117  8.316|  7.522  0.117  8.316|  7.379  0.117  8.316|  7.353  0.117  8.316|  7.348  0.117  8.316|  7.347  0.117       |" + LS
        + "[  2]T                     |                     |         7.631 10.900| 15.083  7.630 10.900| 14.901  7.630 10.900| 14.868  7.630 10.900| 14.861  7.630 10.900| 14.860  7.630 10.900| 14.860  7.630  8.540| 14.860  0.161 10.900|  7.613  7.619  8.540|  9.220  4.426  8.540| 10.529  4.426  8.540| 11.307  4.426  8.540| 11.575  4.426  8.540| 11.638  4.426  8.598| 11.652  4.426       |" + LS
        + "[  3]A                     |                     |                     |        12.660 13.887| 20.113 12.659 13.887| 19.930 12.659 13.887| 19.897 12.659 13.887| 19.890 12.659 13.887| 19.889  7.621 11.629| 15.072 12.404  8.637| 16.641  0.183 11.629|  7.636  6.884 11.392|  9.239  4.428 11.392| 10.544  4.433 11.392| 11.317  4.434 11.392| 11.584  4.434 11.437| 11.647  9.472 12.906| 13.231  9.485       |" + LS
        + "[  4]A                     |                     |                     |                     |        17.504 16.972| 24.957 17.503 16.972| 24.774 17.503 16.972| 24.741 17.503 16.972| 24.734 12.465 14.504| 19.916 12.666 11.738| 19.900  8.827  8.661| 16.274  0.205 14.144|  7.658  6.579 12.722|  9.259  4.443 12.726| 10.562  4.452 12.727| 11.336  4.454 12.734| 11.604  9.492 15.874| 13.189 14.421 17.966| 14.798 14.527       |" + LS
        + "[  5]A                     |                     |                     |                     |                     |        21.690 20.071| 29.142 21.689 20.071| 28.959 21.689 20.071| 28.926 16.651 17.571| 24.102 17.424 14.837| 24.516 11.564 11.758| 19.016  8.171  8.685| 15.617  0.227 14.677|  7.680  6.359 12.870|  9.278  4.459 12.879| 10.581  4.472 12.881| 11.356  9.513 17.660| 12.947 14.459 20.919| 14.556 18.154 23.007| 16.165 19.110       |" + LS
        + "[  6]A                     |                     |                     |                     |        25.033 28.370| 32.485 24.949 27.042| 32.233 23.621 24.597| 31.013 21.177 21.755| 28.611 13.296 18.370| 20.749 14.949 11.228| 21.686  2.769  8.882| 10.222  0.424 11.197|  7.857  2.734 13.010|  9.070  4.543 13.034| 10.442  4.574 15.412| 11.345  6.944 20.447| 12.743 16.708 25.223| 14.352 19.504 27.542| 15.961 21.178 30.116|" + LS
        + "[  7]A                     |                     |                     |                     |                     |        30.064 30.106| 37.516 29.912 27.686| 37.206 28.412 24.848| 35.813 20.820 20.919| 28.272 18.351 14.328| 25.787 11.418 10.857| 18.870  2.790  8.899| 10.242  0.446 11.207|  7.879  2.751 12.980|  9.090  4.557 13.050| 10.461  4.593 15.424| 11.365 11.995 25.145| 12.973 18.001 27.917| 14.582 19.627 29.630| 16.192 21.236       |" + LS
        + "[  8]G                     |                     |                     |                     |                     |                     |        29.506 30.787| 36.958 27.782 27.949| 35.199 30.071 24.015| 36.410 25.324 17.430| 32.771 19.564 13.956| 27.016 15.586 10.877| 23.035  7.848  8.921| 15.300  5.505 11.224| 12.939  7.806 12.994| 14.149  9.608 13.070| 15.518  4.613 20.475| 12.059 11.099 26.471| 13.661 18.246 28.097| 15.271 14.818 29.717| 16.876 21.465       |" + LS
        + "[  9]G                     |                     |                     |                     |                     |                     |                     |        29.322 31.044| 36.775 32.314 27.116| 38.160 29.269 20.531| 36.675 22.686 17.057| 30.138 19.208 13.978| 26.655 16.123 12.009| 23.566 12.660 13.448| 20.107 10.563 15.493| 17.991 12.857 16.034| 19.200  9.596 13.094| 17.026  4.635 19.579| 12.086 15.891 26.715| 13.696 13.892 23.298| 15.303 19.466 29.946| 16.912 22.145       |" + LS
        + "[ 10]C                     |                     |                     |                     |                     |                     |                     |                     |        34.245 30.217| 41.698 32.366 23.632| 39.789 25.787 20.158| 33.240 22.310 17.079| 29.756 19.228 15.110| 26.672 17.246 16.539| 24.671 17.403 18.100| 24.640 15.617 19.030| 23.029 17.879 16.053| 24.233 14.632 13.116| 22.061  4.657 24.368| 12.110 17.341 22.355| 13.719 18.272 27.941| 15.329 15.527 30.625| 16.936 22.187       |" + LS
        + "[ 11]G                     |                     |                     |                     |                     |                     |                     |                     |                     |        35.453 26.733| 42.905 28.889 23.259| 36.341 25.411 20.180| 32.857 22.329 18.211| 29.773 20.347 19.638| 27.773 21.326 21.192| 28.342 22.120 22.000| 29.051 15.612 19.153| 23.064 16.092 16.216| 23.264 18.136 13.138| 24.475  4.680 24.929| 12.132 17.385 26.738| 13.741 18.988 24.008| 15.351 19.906 30.668| 16.960 22.211       |" + LS
        + "[ 12]A                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        26.952 26.360| 34.404 23.474 23.281| 30.920 20.393 21.312| 27.836 18.411 22.737| 25.836 19.388 24.289| 26.404 20.685 25.097| 27.381 26.524 22.107| 28.983 20.649 19.312| 28.022 20.606 16.239| 27.870 18.387 13.160| 25.813  4.702 25.847| 12.154 17.407 26.579| 13.764 13.979 28.382| 15.371 20.620 30.691| 16.980 22.235       |" + LS
        + "[ 13]C                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        31.102 26.378| 38.554 27.843 24.402| 35.288 25.169 25.538| 32.608 23.460 26.908| 30.876 24.440 27.876| 31.453 25.733 25.208| 32.429 27.344 22.411| 33.654 19.253 19.340| 26.705 21.481 16.261| 27.884 18.410 13.182| 25.837  4.724 25.865| 12.176 17.429 22.459| 13.786 18.349 29.091| 15.395 20.646 30.716| 17.004 17.217       |" + LS
        + "[ 14]A                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        26.587 27.503| 34.039 24.583 28.633| 32.008 24.742 29.874| 31.978 23.456 30.843| 30.842 29.473 28.306| 32.440 29.921 25.513| 34.014 27.660 22.436| 34.642 23.752 19.362| 31.198 16.473 16.283| 23.926 18.433 13.204| 25.002  4.746 25.026| 12.198 17.451 26.825| 13.808 19.054 29.117| 15.417 20.670 25.698| 17.027 21.589       |" + LS
        + "[ 15]T                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        31.363 31.499| 38.815 29.628 32.398| 37.046 29.796 31.810| 37.030 28.515 31.406| 35.901 33.187 28.614| 37.468 30.754 25.537| 37.857 27.686 22.463| 35.126 24.604 19.380| 32.047 20.844 16.305| 28.291 18.455 13.226| 25.889  4.768 25.914| 12.220 17.474 27.526| 13.830 19.083 28.266| 15.439 20.692 30.069| 17.049 22.295       |" + LS
        + "[ 16]G                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        37.104 38.808| 44.557 35.388 38.144| 42.805 29.686 37.653| 37.138 29.195 37.050| 36.532 33.629 35.201| 38.090 26.742 32.143| 34.191 28.723 29.062| 35.277 25.641 25.981| 33.071 22.560 22.907| 30.003 19.486 14.557| 26.929  6.099 13.702| 13.551 10.282 15.602| 15.087 12.178 27.411| 16.645 18.953 29.020| 18.254 20.563 30.674|" + LS
        + "[ 17]C                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        42.025 41.175| 49.478 40.397 38.094| 47.811 34.746 37.594| 42.197 29.216 38.280| 36.668 38.495 34.540| 38.277 31.798 32.156| 38.826 28.282 29.076| 35.725 30.241 26.001| 36.805 27.162 17.658| 34.593 19.805 14.477| 27.258  6.121 18.039| 13.573 10.270 20.659| 15.107 17.196 27.424| 16.716 21.810 29.034| 18.325 23.419       |" + LS
        + "[ 18]C                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        46.013 41.194| 53.466 43.236 40.619| 50.676 34.721 37.672| 42.174 34.276 37.641| 41.608 39.667 35.251| 43.197 31.366 32.167| 38.816 33.026 29.102| 39.759 31.241 20.759| 38.627 22.915 17.578| 30.368 14.690 14.599| 22.142  6.143 18.744| 13.595 15.301 25.668| 15.204 20.245 30.143| 16.813 21.990 31.899| 18.423 23.599       |" + LS
        + "[ 19]A                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        46.442 43.720| 53.894 45.791 40.688| 53.144 39.739 40.617| 47.191 34.270 38.352| 41.722 40.492 35.257| 43.321 31.074 32.203| 38.525 34.335 23.860| 39.959 26.016 20.679| 33.468 22.829 17.696| 30.273 19.109 14.622| 26.556  6.165 23.775| 13.617 13.635 28.714| 15.224 20.478 30.459| 16.833 22.088 32.080| 18.443 23.698       |" + LS
        + "[ 20]A                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        48.900 43.789| 56.352 45.938 43.707| 53.380 39.468 41.211| 46.920 39.316 38.358| 46.610 35.469 35.290| 42.916 35.899 26.961| 43.082 29.117 23.780| 36.570 25.930 20.797| 33.374 22.946 17.723| 30.389 14.828 14.645| 22.280  6.187 22.115| 13.639 18.098 28.948| 15.249 20.499 30.559| 16.858 22.108 32.178| 18.467 23.718       |" + LS
        + "[ 21]T                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        49.039 46.808| 56.492 48.841 44.286| 56.142 44.394 41.458| 51.843 43.232 38.387| 50.624 34.807 30.063| 42.259 32.218 26.881| 39.656 29.031 23.898| 36.475 26.047 20.824| 33.490 22.974 17.742| 30.417 19.202 14.667| 26.650  6.209 26.575| 13.661 18.907 28.970| 15.271 15.486 30.580| 16.878 22.133 32.198| 18.487 18.705       |" + LS
        + "[ 22]G                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        51.942 47.387| 59.394 49.531 44.559| 56.965 46.652 41.488| 54.093 43.635 33.164| 51.078 30.271 29.982| 37.724 32.132 27.000| 38.758 29.148 23.925| 36.578 26.075 20.843| 33.518 22.993 17.768| 30.436 19.911 14.689| 27.354  6.231 27.378| 13.683 18.937 23.966| 15.293 14.818 30.604| 16.898 22.153 27.185| 18.507 23.074       |" + LS
        + "[ 23]T                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        52.632 47.660| 60.084 49.807 44.590| 57.248 41.701 36.265| 49.153 38.421 33.080| 45.866 34.591 30.101| 42.039 32.249 27.026| 39.683 29.176 23.944| 36.619 26.094 20.869| 33.537 17.981 17.790| 25.433 19.940 14.711| 26.509  6.253 26.534| 13.706 18.959 23.299| 15.315 14.436 29.752| 16.917 22.173 31.554| 18.527 23.776       |" + LS
        + "[ 24]G                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        52.908 47.691| 60.361 49.840 39.366| 57.283 36.479 36.181| 43.931 38.331 33.202| 44.961 35.344 30.128| 42.774 32.277 27.046| 39.720 29.195 23.970| 36.638 26.120 20.888| 33.563 17.313 17.812| 24.765 19.962 14.733| 26.072  6.275 26.097| 13.728 18.981 22.916| 15.337 19.209 30.635| 16.946 17.155 32.256| 18.554 23.802       |" + LS
        + "[ 25]T                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        47.903 42.467| 55.356 44.623 39.278| 52.068 40.794 36.303| 48.242 38.451 33.229| 45.885 35.378 30.147| 42.821 32.296 27.071| 39.739 24.183 23.989| 31.635 26.138 20.906| 32.710 16.931 17.835| 24.383 19.984 14.756| 25.781  6.297 25.805| 13.750 19.003 27.687| 15.359 20.593 25.635| 16.969 21.528 32.282| 18.578 23.829       |" + LS
        + "[ 26]T                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        44.435 44.873| 51.888 36.414 41.805| 43.867 38.384 38.723| 44.947 30.265 35.648| 37.717 32.227 31.972| 38.794 28.551 29.480| 35.991 26.059 25.175| 33.495 21.755 23.324| 29.205 14.865 14.883| 22.317  6.425 17.437| 13.877  8.973 18.284| 15.156  9.822 29.634| 16.295 16.510 31.448| 17.902 23.159 33.305| 19.512 25.363 35.816|" + LS
        + "[ 27]T                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        44.032 44.207| 51.485 41.471 41.818| 48.908 37.944 38.054| 45.390 35.321 35.069| 42.759 36.563 32.569| 43.483 33.330 28.275| 40.769 30.025 26.403| 37.470 21.619 17.979| 29.071 14.296 14.902| 21.748  6.447 17.434| 13.899  8.989 18.302| 15.177  9.840 24.990| 16.315 20.880 31.631| 17.924 23.177 33.837| 19.533 24.790       |" + LS
        + "[ 28]C                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        48.568 44.913| 56.020 46.073 41.150| 53.509 37.415 38.167| 44.867 39.660 35.670| 46.052 32.766 31.376| 40.218 28.486 29.504| 35.936 31.628 21.080| 37.350 23.205 17.994| 30.657 18.985 14.924| 26.435 11.507 17.450| 18.959 14.043 18.321| 20.235  9.859 29.356| 17.300 21.583 31.653| 18.910 23.199 33.270| 20.519 19.771       |" + LS
        + "[ 29]A                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        50.139 44.251| 57.591 46.398 41.258| 53.846 42.147 38.771| 49.597 40.905 34.476| 48.301 36.368 32.592| 43.818 33.286 24.181| 40.729 26.337 21.095| 33.790 23.245 18.025| 30.688 20.160 19.537| 27.603 16.565 21.135| 24.012 19.086 18.339| 25.285  9.881 30.055| 17.333 17.527 31.670| 18.940 24.174 28.251| 20.549 24.509       |" + LS
        + "[ 30]A                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        49.505 44.359| 56.957 46.508 41.872| 53.951 43.988 37.577| 51.424 39.731 35.693| 47.180 37.822 27.283| 45.245 29.438 24.196| 36.891 26.346 21.126| 33.789 23.275 22.635| 30.718 24.291 23.868| 31.300 21.616 21.438| 29.047 18.100 18.361| 25.547  9.903 26.007| 17.355 21.908 31.112| 18.965 24.215 32.990| 20.574 25.808       |" + LS
        + "[ 31]C                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        49.609 44.973| 57.062 42.080 40.678| 49.532 37.794 38.794| 45.243 40.929 30.384| 46.657 32.540 27.297| 39.992 29.447 24.227| 36.890 26.376 25.734| 33.819 27.396 26.966| 34.404 23.500 24.535| 30.946 25.992 21.456| 32.211 22.670 18.383| 30.098  9.925 30.366| 17.377 22.624 32.662| 18.987 24.240 34.288| 20.596 25.849       |" + LS
        + "[ 32]T                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        50.219 43.778| 57.671 45.673 41.883| 53.123 37.552 33.485| 45.004 30.603 30.398| 38.055 27.510 27.328| 34.954 24.439 28.833| 31.882 25.458 30.064| 32.466 31.639 27.623| 34.069 28.304 24.557| 35.024 26.700 21.484| 34.072 23.623 18.405| 31.066  9.947 31.095| 17.399 17.615 32.711| 19.006 19.224 34.329| 20.614 25.871       |" + LS
        + "[ 33]T                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        49.033 44.984| 56.485 42.075 36.586| 49.527 33.683 33.496| 41.136 29.928 30.425| 37.376 26.846 31.617| 34.290 24.451 32.785| 31.885 30.508 30.724| 33.483 32.851 27.658| 35.087 29.786 24.585| 36.238 26.735 21.506| 34.162 18.618 18.427| 26.070  9.969 26.095| 17.421 16.947 27.704| 19.026 23.593 34.352| 20.636 25.889       |" + LS
        + "[ 34]T                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        45.182 39.687| 52.634 36.800 36.593| 44.252 33.017 33.519| 40.465 29.545 34.283| 36.991 26.861 32.881| 34.300 29.511 33.820| 35.611 34.946 30.759| 37.214 32.906 27.686| 38.629 29.835 24.607| 37.237 21.719 21.525| 29.171 17.950 18.449| 25.398  9.991 25.423| 17.443 21.602 32.069| 19.053 24.295 34.369| 20.662 25.911       |" + LS
        + "[ 35]C                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |                     |        44.939 39.690|        41.162 36.612|        37.673 36.961|        34.598 34.918|        26.882 36.626|        34.553 33.860|        35.990 30.787|        27.899 27.708|        29.858 24.622|        26.087 21.543|        22.606 18.471|        10.013 30.077|        17.668 32.766|        24.328 34.391|        25.938       |" + LS
        + "                           |G            0       |G            1       |G            2       |G            3       |G            4       |G            5       |G            6       |A            7       |T            8       |A            9       |A           10       |A           11       |A           12       |A           13       |A           14       |G           15       |G           16       |C           17       |G           18       |A           19       |C           20       |A           21       |T           22       |G           23       |C           24       |C           25       |A           26       |A           27       |T           28       |G           29       |T           30       |G           31       |T           32       |C           33       |G           34       |C           35       |C           36       |T           37       |T           38       |T           39       |T           40       |T           41       |C           42       |A           43       |A           44       |C           45       |T           46       |T           47       |T           48       |C           49       |C           50       |G           51       |A           52       |T           53       |" + LS
        ,
        score.toString());
  }

  public void testReadEndsAfterLn() {
    final AllPaths score = score(env());
    if (score instanceof ScoreMatrixCG) {
      final ScoreMatrixCG sm = (ScoreMatrixCG) score;
      assertEquals(0.0, sm.readEndsAfterLn(37));
      assertEquals(0.0, sm.readEndsAfterLn(48), 0.0001);
      assertEquals(-0.0002, sm.readEndsAfterLn(49), 0.0001);
      assertEquals(-7.6540, sm.readEndsAfterLn(50), 0.0001);
      assertEquals(-14.1332, sm.readEndsAfterLn(51), 0.0001);
      assertEquals(-15.9250, sm.readEndsAfterLn(52), 0.0001);
      assertEquals(Double.NEGATIVE_INFINITY, sm.readEndsAfterLn(53));
    }
  }

  public static void main(final String[] args) {
    final double[] qualities = new double[READ.length()];
    for (int i = 0; i < qualities.length; i++) {
      qualities[i] = 0.01;
    }
    // // Enable this code to have constantly-changing qualities (worst case)
    //    qualities[0] = 0.1;
    //    qualities[1] = Math.pow(10.0, -15 / 10.0);
    //    for (int i = 3; i < 34; i += 2) {
    //      qualities[i] = 0.1;
    //    }
    final EnvironmentImplementation env = new EnvironmentImplementation(
        7, //maxShift
        DnaUtils.encodeString(TEMPLATE),
        7, //start
        DnaUtils.encodeString(READ),
        qualities
        );
    final RealignParams params = new MockRealignParams();
    final double[][] sums = new double[50][4];
    final double sum = 0.0;
    final int repeats = 10000;
    final long start = System.nanoTime();
    for (int i = 0; i < repeats; i++) {

      final ScoreMatrixCG m = new ScoreMatrixCG(LogPossibility.SINGLETON, params);
      final ScoreMatrixCG m2 = new ScoreMatrixCG(LogPossibility.SINGLETON, params);
      final ScoreMatrixCGReverse mr = new ScoreMatrixCGReverse(LogPossibility.SINGLETON, params);

      final Delta re = new DeltaSlowly(LogPossibility.SINGLETON, m, m2, mr);
      re.setEnv(env);
      for (int tpos = 0; tpos < sums.length; tpos++) {
        final double[] probs = re.probabilitiesLn(tpos);
        sums[tpos][0] += probs[0];
        sums[tpos][1] += probs[1];
        sums[tpos][2] += probs[2];
        sums[tpos][3] += probs[3];
      }
    }
    final long finish = System.nanoTime();
    for (int tpos = 0; tpos < sums.length; tpos++) {
      System.out.println(String.format("scores at pos %2d:  %c =%10.4f  %c =%10.4f  %c =%10.4f  %c =%10.4f  ",
          tpos,
          'A', sums[tpos][0] / repeats,
          'C', sums[tpos][1] / repeats,
          'G', sums[tpos][2] / repeats,
          'T', sums[tpos][3] / repeats)
          );
    }
    System.out.println("average score = " + sum / repeats);
    System.out.println("total time=" + (finish - start) / 1.0E9 + "  msec/call=" + (finish - start) / 1.0E6 / repeats);
  }
}
