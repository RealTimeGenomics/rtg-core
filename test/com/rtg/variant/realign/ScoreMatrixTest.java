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

import java.io.IOException;
import java.util.Arrays;

import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 */
public class ScoreMatrixTest extends TestCase {

  static final String HOMOPOLYMER_EXP = ""
      + "ScoreMatrix                |G            0       |G            1       |A            2       |C            3       |G            4       |G            5       |G            6       |C            7       |G            8       |T            9       |T           10       |T           11       |N           12       |" + LS
      + "[  0]          0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|" + LS
      + "[  1]A                     |         3.327  8.316| 10.780  3.327  8.316| 10.598  0.117  8.316|  7.560  3.327  8.316|  8.987  3.327       |" + LS
      + "[  2]C                     |                     |         7.631 10.900| 15.083  7.630  8.540| 14.901  0.161 10.900|  7.613  7.619 11.808|  9.220  7.634       |" + LS
      + "[  3]G                     |                     |                     |        12.660 11.629| 20.113 12.404  8.637| 19.713  0.183 14.643|  7.636  7.044 16.114|  9.240  7.504       |" + LS
      + "[  4]G                     |                     |                     |                     |        16.526 11.738| 23.978  8.828  8.663| 16.281  0.205 15.499|  7.658  6.699 15.984|  9.260 12.432       |" + LS
      + "[  5]C                     |                     |                     |                     |                     |        16.985 11.761| 24.437 13.211  8.686| 20.659  5.265 15.159| 12.718  6.453 20.912| 13.401 14.490       |" + LS
      + "[  6]G                     |                     |                     |                     |                     |                     |        11.973 11.787| 19.426  8.892 13.735| 16.335 10.325 14.933| 17.164  6.474 22.971| 13.918 18.332       |" + LS
      + "[  7]T                     |                     |                     |                     |                     |                     |                     |        16.345 16.375| 23.798 13.945 17.654| 21.380 15.375 14.954| 22.212  6.496 26.813| 13.948 14.147       |" + LS
      + "[  8]T                     |                     |                     |                     |                     |                     |                     |                     |        20.819 20.583|        18.985 18.052|        14.585 14.976|         6.518 22.627|        13.487       |" + LS
      + "                           |G            0       |G            1       |A            2       |C            3       |G            4       |G            5       |G            6       |C            7       |G            8       |T            9       |T           10       |T           11       |N           12       |" + LS
      ;
  //Mark says I have to put my name here so he doesnt get blamed for the next bit of code
  //JC

  /** Used for testing. */
  public static class LogTest extends ScoreMatrixTest {
    public LogTest() { }
    @Override
    protected AllPaths score(final Environment env) {
      return score(LogPossibility.SINGLETON, env);
    }
  }

  /** Used for testing. */
  public static class LogApproximateTest extends ScoreMatrixTest {
    public LogApproximateTest() { }
    @Override
    protected AllPaths score(final Environment env) {
      return score(LogApproximatePossibility.SINGLETON, env);
    }
    //the results are slightly inaccurate in the 5th digit - but the other two accurate versions agree so let it be
    @Override
    public void test1() { }

    @Override
    public void test2() { }
  }

  /** Used for testing. */
  public static class SimpleTest extends ScoreMatrixTest {
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
   * Realign parameters that match the scorematrixtestCG.xls spreadsheet.
   */
  public static class MockRealignParams extends IntegralAbstract implements RealignParams {
    @Override
    public double matchLn() {
      return Math.log(0.98939);
    }
    @Override
    public double misMatchLn() {
      return Math.log(0.00920);
    }
    @Override
    public double insertOpenLn() {
      return Math.log(0.00058);
    }
    @Override
    public double insertExtendLn() {
      return Math.log(0.2);
    }
    @Override
    public double deleteOpenLn() {
      return Math.log(0.00083);
    }
    @Override
    public double deleteExtendLn() {
      return Math.log(0.18);
    }
    @Override
    public boolean completeGenomics() {
      return false;
    }
    @Override
    public int gapEnd(final int gap) {
      throw new UnsupportedOperationException();
    }
    @Override
    public double gapFreqLn(final int gap, final int width) {
      throw new UnsupportedOperationException();
    }
    @Override
    public int gapStart(final int gap) {
      throw new UnsupportedOperationException();
    }
    @Override
    public void toString(final StringBuilder sb) {
      // do nothing
    }
    @Override
    public boolean integrity() {
      return true;
    }
    @Override
    public double[][] gapDistributionPoss(PossibilityArithmetic arith) {
      throw new UnsupportedOperationException();
    }
  }

  static final String READ = "ACGCGACG";
  static final String TEMPLATE = "GGACGTACGTTT";

  static final String READ_H = "ACGGCGTT";
  static final String TEMPLATE_H = "GGACGGGCGTTT";

  protected Environment env() {
    return env(READ, TEMPLATE, 2);
  }

  protected Environment env(String read, String template, int maxShift) {
    //see scorematrixtest.xls for details
    final int start = 2;
    final byte[] readBytes = DnaUtils.encodeString(read);
    final double[] quality = new double[readBytes.length];
    Arrays.fill(quality, 0.01);
    quality[0] = 0.1;
    quality[1] = Math.pow(10.0, -15 / 10.0);
    final byte[] templateBytes = DnaUtils.encodeString(template);
    final Environment env = new EnvironmentImplementation(
        maxShift,
        templateBytes,
        start,
        readBytes,
        quality
        );
    Exam.integrity(env);
    return env;
  }

  protected AllPaths score(final Environment env) {
    return score(LogPossibility.SINGLETON, env);
  }

  protected AllPaths score(final PossibilityArithmetic arith, final Environment env) {
    final AllPaths score = new ScoreMatrix(arith, new MockRealignParams());
    score.setEnv(env);
    Exam.globalIntegrity(score);
    return score;
  }

  private static final String READ_EXP = ""
      + "ScoreMatrix                |G            0       |G            1       |A            2       |C            3       |G            4       |T            5       |A            6       |C            7       |G            8       |T            9       |T           10       |T           11       |N           12       |" + LS
      + "[  0]          0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|" + LS
      + "[  1]A                     |         3.327  8.316| 10.780  3.327  8.316| 10.598  0.117  8.316|  7.560  3.327  8.316|  8.987  3.327       |" + LS
      + "[  2]C                     |                     |         7.631 10.900| 15.083  7.630  8.540| 14.901  0.161 10.900|  7.613  7.619 11.808|  9.220  7.634       |" + LS
      + "[  3]G                     |                     |                     |        12.660 11.629| 20.113 12.404  8.637| 19.713  0.183 14.643|  7.636 12.082 16.114|  9.245 12.542       |" + LS
      + "[  4]C                     |                     |                     |                     |        11.488 11.738| 18.940 13.866  8.663| 20.169  5.243 18.984| 12.696 12.903 21.022| 14.303  9.444       |" + LS
      + "[  5]G                     |                     |                     |                     |                     |        11.015 11.765| 18.468 13.914 13.723| 19.834 10.303 21.321| 17.731 17.277 17.924| 19.336  9.460       |" + LS
      + "[  6]A                     |                     |                     |                     |                     |                     |        15.747 16.821| 23.200 13.242 18.780| 20.678 15.363 21.016| 21.824 21.675 17.940| 23.430 14.520       |" + LS
      + "[  7]C                     |                     |                     |                     |                     |                     |                     |        20.560 21.105| 28.012 13.261 23.278| 20.713 20.416 21.041| 22.319 23.149 23.000| 23.927 19.580       |" + LS
      + "[  8]G                     |                     |                     |                     |                     |                     |                     |                     |        25.231 21.731|        13.283 24.134|        24.766 26.097|        26.876 28.060|        24.629       |" + LS
      + "                           |G            0       |G            1       |A            2       |C            3       |G            4       |T            5       |A            6       |C            7       |G            8       |T            9       |T           10       |T           11       |N           12       |" + LS
      ;

  //typical case
  public void test1() throws IOException {
    final AllPaths score = score(env());
    assertEquals(READ_EXP, score.toString());
  }

  //homopolymer sequence for checking against tests in HomopolymerMatrixTest
  public void test2() throws IOException {
    final AllPaths score = score(env(READ_H, TEMPLATE_H, 2));
    final String exp = HOMOPOLYMER_EXP
        ;
    assertEquals(exp, score.toString());
  }

  /*
  public void testTotalScoreNs() throws IOException {
    for (String read : new String[]{"GGGA", "ACCG", "CG"}) {
      for (String template : new String[] {"NNNN", "GGGA", "ACCG", "CG"}) {
        final Environment env = env(read, template, 2);
        final AllPaths score = score(env);
        System.err.print(read + "->" + template + " " + score.totalScore() + "   ");
      }
      System.err.println();
    }
  }
  */

  public void testTotalScoreLn() throws IOException {
    final Environment env = env();
    final AllPaths score = score(env);
    final double expected = -13.2824;
    assertEquals(expected, score.totalScoreLn(), 0.0001);
    assertEquals(expected, Math.log(score.totalScore()), 0.0001);
    assertEquals(expected, score.arithmetic().poss2Ln(score.total()), 0.0001);
    if (score instanceof AbstractAllPaths) {
      assertEquals(-env.maxShift(), ((AbstractAllPaths) score).rowOffset(1));
    }
  }

  public void testReadEndsAfterLn() throws IOException {
    final AllPaths score = score(env());
    if (score instanceof ScoreMatrix) {
      final ScoreMatrix sm = (ScoreMatrix) score;
      assertEquals(0.0, sm.readEndsAfterLn(-1));
      assertEquals(0.0, sm.readEndsAfterLn(7));
      assertEquals(-0.0002, sm.readEndsAfterLn(8), 0.0001);
      assertEquals(-10.5400, sm.readEndsAfterLn(9), 0.0001);
      assertEquals(-11.2175, sm.readEndsAfterLn(10), 0.0001);
      assertEquals(-11.3468, sm.readEndsAfterLn(11), 0.0001);
      assertEquals(Double.NEGATIVE_INFINITY, sm.readEndsAfterLn(12));
      assertFalse(sm.underflow());
    }
  }

  private static final String SMALL_READ = "ACGCA";

  private static final String SMALL_EXP = ""
      + "ScoreMatrix                |G            0       |G            1       |A            2       |C            3       |G            4       |T            5       |A            6       |C            7       |G            8       |T            9       |" + LS
      + "[  0]          0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|" + LS
      + "[  1]A                     |         3.327  8.316| 10.780  3.327  8.316| 10.598  0.117  8.316|  7.560  3.327  8.316|  8.987  3.327       |" + LS
      + "[  2]C                     |                     |         7.631 10.900| 15.083  7.630  8.540| 14.901  0.161 10.900|  7.613  7.619 11.808|  9.220  7.634       |" + LS
      + "[  3]G                     |                     |                     |        12.660 11.629| 20.113 12.404  8.637| 19.713  0.183 14.643|  7.636 12.082 16.114|  9.245 12.542       |" + LS
      + "[  4]C                     |                     |                     |                     |        11.488 11.738| 18.940 13.866  8.663| 20.169  5.243 18.984| 12.696 12.903 21.022| 14.303  9.444       |" + LS
      + "[  5]A                     |                     |                     |                     |                     |        16.053 11.765|        13.914 13.723|         5.265 21.321|        17.277 17.924|        14.498       |" + LS
      + "                           |G            0       |G            1       |A            2       |C            3       |G            4       |T            5       |A            6       |C            7       |G            8       |T            9       |" + LS;


  public void testSmallMatrix() {
    final AllPaths score = new ScoreMatrix(LogPossibility.SINGLETON, new MockRealignParams());
    score.setEnv(env(SMALL_READ, TEMPLATE, 2));
    Exam.globalIntegrity(score);

    assertEquals(-5.263277, score.totalScoreLn(), 0.000001);
    assertEquals(SMALL_EXP, score.toString());
  }

  public void testMatrixReuse() {

    final int[] envcalls = new int[1];
    final int[] resizecalls = new int[1];

    final AbstractAllPaths score = new ScoreMatrix(LogPossibility.SINGLETON, new MockRealignParams()) {

      @Override
      public void setEnv(final Environment env) {
        super.setEnv(env);
        envcalls[0]++;
      }

      @Override
      public void resizeMatrix(int length, int width) {
        super.resizeMatrix(length, width);
        resizecalls[0]++;
      }
    };

    assertEquals(0, envcalls[0]);
    assertEquals(0, resizecalls[0]);

    score.setEnv(env(READ, TEMPLATE, 4));
    score.globalIntegrity();
    final String exp = ""
        + "ScoreMatrix                |N           -2       |N           -1       |G            0       |G            1       |A            2       |C            3       |G            4       |T            5       |A            6       |C            7       |G            8       |T            9       |T           10       |T           11       |N           12       |N           13       |N           14       |" + LS
        + "[  0]          0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|" + LS
        + "[  1]A                     |         1.388  8.316|  8.840  1.388  8.316|  8.658  3.327  8.316|  9.798  3.327  8.316| 10.352  0.117  8.316|  7.557  3.327  8.316|  8.985  3.327  8.316|  9.990  3.327  8.316| 10.415  0.117       |" + LS
        + "[  2]C                     |                     |         2.775  9.675| 10.227  5.696 10.900| 11.598  7.627 10.900| 13.065  7.630  8.540| 14.164  0.161 10.900|  7.613  7.619 10.900|  9.220  7.628 10.900| 10.815  7.630  8.598| 12.357  0.161       |" + LS
        + "[  3]G                     |                     |                     |         2.796 13.392| 10.249  5.705 13.886| 11.617 12.642 11.629| 13.225 12.402  8.637| 14.828  0.183 13.885|  7.636 12.071 13.886|  9.245 12.511 11.687| 10.854 12.393  8.642| 12.463  0.183       |" + LS
        + "[  4]C                     |                     |                     |                     |         7.856 14.126| 15.309 10.756 14.729| 16.675 10.960 11.738| 17.653 13.856  8.663| 19.141  5.243 16.960| 12.695 12.901 14.786| 14.302  9.360 11.743| 15.571 13.772  8.664| 17.163  5.243       |" + LS
        + "[  5]G                     |                     |                     |                     |                     |        12.914 17.611| 20.367 15.793 14.829| 21.729 10.660 11.764| 18.108 13.914 13.722| 19.541 10.303 17.857| 17.723 17.215 14.795| 19.327  9.304 11.765| 16.741 13.913 13.724| 18.303 10.303       |" + LS
        + "[  6]A                     |                     |                     |                     |                     |                     |        17.967 17.928| 25.420 19.702 14.852| 26.396 15.480 16.819| 22.926 13.239 18.676| 20.670 15.363 17.896| 21.819 19.903 14.813| 23.409 14.296 16.821| 21.712 18.260 18.784| 23.233 15.363       |" + LS
        + "[  7]C                     |                     |                     |                     |                     |                     |                     |        22.410 17.953| 29.862 20.099 19.903| 27.532 20.345 21.055| 27.566 13.258 20.940| 20.710 20.356 17.914| 22.315 20.062 19.866| 23.897 19.293 21.877| 25.252 22.907 23.844| 26.832 16.751       |" + LS
        + "[  8]G                     |                     |                     |                     |                     |                     |                     |                     |        23.196 23.000|        24.466 24.146|        25.065 21.643|        13.279 21.015|        23.018 22.964|        24.388 24.919|        24.285 26.933|        23.959 25.231|        18.138       |" + LS
        + "                           |N           -2       |N           -1       |G            0       |G            1       |A            2       |C            3       |G            4       |T            5       |A            6       |C            7       |G            8       |T            9       |T           10       |T           11       |N           12       |N           13       |N           14       |" + LS
        ;
    assertEquals(exp, score.toString());
    assertEquals(-13.2705, score.totalScoreLn(), 0.0001);

    assertEquals(1, envcalls[0]);
    assertEquals(1, resizecalls[0]);

    score.setEnv(env(SMALL_READ, TEMPLATE, 2));
    score.globalIntegrity();
    assertEquals(-5.263277, score.totalScoreLn(), 0.000001);
    assertEquals(SMALL_EXP, score.toString());

    assertEquals(2, envcalls[0]);
    assertEquals(1, resizecalls[0]);
  }

  /**
   * For speed testing with various arithmetic implementations.
   * You need to make this class non-abstract to run this.
   *
   * @param args none required.
   */
  public static void main(String[] args) {
    final String read = "ACGTACGTAC GTACGTACGT ACGTACGTAC GTACG".replaceAll(" ", "");
    final String tmpl = "AACCGGTTAACCGGTTTGCAGACTGACTGACGAGTACATGCTGCAGCGTAGCG";
    final int start = 10;
    final int maxShift = 7;
    final int repeat = 1000000;
    final double[] quality = new double[read.length()];
    Arrays.fill(quality, 0.01);
    final Environment env = new EnvironmentImplementation(
        maxShift,
        DNA.stringDNAtoByte(tmpl),
        start,
        DNA.stringDNAtoByte(read),
        quality
        );
    Exam.integrity(env);
    final long time0 = System.nanoTime();
    ScoreMatrix sm1 = null;
    for (int i = 0; i < repeat; i++) {
      sm1 = new ScoreMatrix(LogPossibility.SINGLETON, new MockRealignParams());
      sm1.setEnv(env);
    }
    final long time1 = System.nanoTime();
    ScoreMatrix sm2 = null;
    for (int i = 0; i < repeat; i++) {
      sm2 = new ScoreMatrix(LogApproximatePossibility.SINGLETON, new MockRealignParams());
      sm2.setEnv(env);
    }
    final long time2 = System.nanoTime();
    ScoreMatrix sm3 = null;
    for (int i = 0; i < repeat; i++) {
      sm3 = new ScoreMatrix(SimplePossibility.SINGLETON, new MockRealignParams());
      sm3.setEnv(env);
    }
    final long time3 = System.nanoTime();
    final double usec1 = (time1 - time0) / (1000.0 * repeat);
    final double usec2 = (time2 - time1) / (1000.0 * repeat);
    final double usec3 = (time3 - time2) / (1000.0 * repeat);
    System.out.println("log       scoreLn=" + sm1.totalScoreLn() + "  time=" + usec1 + "usec");
    System.out.println("logApprox scoreLn=" + sm2.totalScoreLn() + "  time=" + usec2 + "usec");
    System.out.println("double    scoreLn=" + sm3.totalScoreLn() + "  time=" + usec3 + "usec");
    System.out.println("ratios = " + Utils.realFormat(usec1 / usec3, 1)
        + " : " + Utils.realFormat(usec2 / usec3, 1) + " : 1.0");
  }
}
