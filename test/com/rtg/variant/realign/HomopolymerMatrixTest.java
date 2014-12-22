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
import java.io.Reader;
import java.io.StringReader;
import java.util.Arrays;

import com.rtg.mode.DnaUtils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

import junit.framework.TestCase;

/**
 */
public class HomopolymerMatrixTest extends TestCase {

  private static final String IONT = ""
      + "A+T 3 0 10 10 10" + LS
      + "A+T 4 0 0 10 10 10" + LS
      + "C+G 3 0 42 42 0" + LS
      + "C+G 4 0 0 0 42 42" + LS
      ;
  protected AllPaths score(final EnvironmentHomopolymer env) throws IOException {
    return score(LogPossibility.SINGLETON, env, IONT);
  }

  protected AllPaths score(final PossibilityArithmetic arith, final EnvironmentHomopolymer env, final String calibration) throws IOException {
    final Reader in = new StringReader(calibration.replaceAll(" ", "\t"));
    final HomoPolymerParams params = new HomoPolymerParams(arith, 2, 3, in);
    //System.err.println(params);
    //System.err.println(params.transition(2, 3, 2));
    final AllPaths score = new HomopolymerMatrix(arith, new MockRealignParams(), params);
    score.setEnv(env);
    Exam.globalIntegrity(score);
    return score;
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

  protected EnvironmentHomopolymer env(final String read, final String template) {
    //see scorematrixtest.xls for details
    final int start = 2;
    final int maxShift = 2;
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
    return new EnvironmentHomopolymer(env);
  }


  public void test1() throws IOException {
    final EnvironmentHomopolymer env = env(ScoreMatrixTest.READ, ScoreMatrixTest.TEMPLATE);
    //System.err.println(env);
    final AllPaths score = score(env);
    final String exp = ""
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
    assertEquals(exp, score.toString());
  }

  //homopolymer reads that invoke short circuit calculations in HomopolymerMatrix
  public void test2a() throws IOException {
    final EnvironmentHomopolymer env = env(ScoreMatrixTest.READ_H, ScoreMatrixTest.TEMPLATE_H);
    //System.err.println(env);
    final AllPaths score = score(LogPossibility.SINGLETON, env, IONT);
    final String exp = ""
        + "ScoreMatrix                |G            0       |G            1       |A            2       |C            3       |G            4       |G            5       |G            6       |C            7       |G            8       |T            9       |T           10       |T           11       |N           12       |" + LS
        + "[  0]          0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|         0.001  7.094|" + LS
        + "[  1]A                     |         3.327  8.316| 10.780  3.327  8.316| 10.598  0.117  8.316|  7.560  3.327  8.316|  8.987  3.327       |" + LS
        + "[  2]C                     |                     |         7.631 10.900| 15.083  7.630  8.540| 14.901  0.161 10.900|  7.613  7.619 11.808|  9.220  7.634       |" + LS
        + "[  3]G                     |                     |                     |        12.660 11.629| 20.113 12.404  8.637| 19.713  0.876 14.643|  8.329  7.044 16.114|  9.928  7.504       |" + LS
        + "[  4]G                     |                     |                     |                     |        16.526 11.738| 23.978  8.828  9.356| 16.281  0.898 15.499|  8.351  3.960 15.984|  9.750 12.495       |" + LS
        + "[  5]C                     |                     |                     |                     |                     |        16.985 12.450| 24.437 13.493  9.379| 20.940  5.958 12.439| 13.411  3.972 20.976| 11.398 14.955       |" + LS
        + "[  6]G                     |                     |                     |                     |                     |                     |        12.656 12.480| 20.108  9.578 14.152| 17.021 11.017 12.453| 17.854  3.994 23.435| 11.446 16.644       |" + LS
        + "[  7]T                     |                     |                     |                     |                     |                     |                     |        17.033 16.884| 24.486 14.630 15.535| 22.064 15.897 12.475| 22.805  4.422 25.125| 11.874 11.683       |" + LS
        + "[  8]T                     |                     |                     |                     |                     |                     |                     |                     |        21.423 18.624|        19.403 15.576|        12.655 12.902|         4.444 20.164|        10.197       |" + LS
        + "                           |G            0       |G            1       |A            2       |C            3       |G            4       |G            5       |G            6       |C            7       |G            8       |T            9       |T           10       |T           11       |N           12       |" + LS
        ;
    assertEquals(exp, score.toString());
  }

  //homopolymer reads that invoke short circuit calculations in HomopolymerMatrix but with calibration table empty so should be same as ScoreMatrix results
  public void test2b() throws IOException {
    final EnvironmentHomopolymer env = env(ScoreMatrixTest.READ_H, ScoreMatrixTest.TEMPLATE_H);
    //System.err.println(env);
    final AllPaths score = score(LogPossibility.SINGLETON, env, "");
    assertEquals(ScoreMatrixTest.HOMOPOLYMER_EXP, score.toString());
  }

  public void testTotalScoreLn() throws IOException {
    final EnvironmentHomopolymer env = env(ScoreMatrixTest.READ_H, ScoreMatrixTest.TEMPLATE_H);
    final AllPaths score = score(LogPossibility.SINGLETON, env, IONT);
    //regression tests
    assertEquals(-4.44028, score.totalScoreLn(), 0.0001);
    assertEquals(-4.44028, Math.log(score.totalScore()), 0.0001);
    if (score instanceof AbstractAllPaths) {
      assertEquals(-env.maxShift(), ((AbstractAllPaths) score).rowOffset(1));
    }
  }

  public void testReadEndsAfterLn() throws IOException {
    final AllPaths score = score(LogPossibility.SINGLETON, env(ScoreMatrixTest.READ_H, ScoreMatrixTest.TEMPLATE_H), IONT);
    if (score instanceof ScoreMatrix) {
      final ScoreMatrix sm = (ScoreMatrix) score;
      //regression tests
      assertEquals(0.0, sm.readEndsAfterLn(-1));
      assertEquals(0.0, sm.readEndsAfterLn(7));
      assertEquals(-7.3387e-7, sm.readEndsAfterLn(8), 0.0001e-7);
      assertEquals(-1.5639e-5, sm.readEndsAfterLn(9), 0.0001e-5);
      assertEquals(-4.9783e-4, sm.readEndsAfterLn(10), 0.0001e-4);
      assertEquals(-5.7572, sm.readEndsAfterLn(11), 0.0001);
      assertEquals(Double.NEGATIVE_INFINITY, sm.readEndsAfterLn(12));
      assertFalse(sm.underflow());
    }
  }

}
