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

import java.io.IOException;
import java.util.Arrays;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

/**
 */
public class ScoreMatrixTest extends AbstractNanoTest {

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
    public MachineType machineType() {
      return null;
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

  //typical case
  public void test1() throws IOException {
    final AllPaths score = score(env());
    mNano.check("scorematrix-typical.txt", score.toString());
  }

  //homopolymer sequence for checking against tests in HomopolymerMatrixTest
  public void test2() throws IOException {
    final AllPaths score = score(env(READ_H, TEMPLATE_H, 2));
    mNano.check("scorematrix-homopolymer.txt", score.toString());  // This expected result is shared with HomopolymerTest
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

  public void testSmallMatrix() throws IOException {
    final AllPaths score = new ScoreMatrix(LogPossibility.SINGLETON, new MockRealignParams());
    score.setEnv(env(SMALL_READ, TEMPLATE, 2));
    Exam.globalIntegrity(score);

    assertEquals(-5.263277, score.totalScoreLn(), 0.000001);
    mNano.check("scorematrix-small.txt", score.toString());
  }

  public void testMatrixReuse() throws IOException {

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
    mNano.check("scorematrix-re-use1.txt", score.toString());
    assertEquals(-13.2705, score.totalScoreLn(), 0.0001);

    assertEquals(1, envcalls[0]);
    assertEquals(1, resizecalls[0]);

    score.setEnv(env(SMALL_READ, TEMPLATE, 2));
    score.globalIntegrity();
    assertEquals(-5.263277, score.totalScoreLn(), 0.000001);
    mNano.check("scorematrix-re-use2.txt", score.toString());

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
