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

import java.io.IOException;
import java.util.Arrays;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.launcher.CommonFlags;
import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.util.MathUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CFlags;
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

  private static final String TEMPLATE_FLAG = "template";
  private static final String READ_FLAG = "read";
  private static final String START_FLAG = "start";
  private static final String MAXSHIFT_FLAG = "maxshift";
  private static final String QUALITY_FLAG = "quality";
  private static final String PPM_FLAG = "ppm";
  private static final String ASCII_FLAG = "ascii";
  private static final String ITERATE_FLAG = "iterate";

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

  public void testTotalScoreLn() {
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

  public void testReadEndsAfterLn() {
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
    final CFlags f = new CFlags();
    f.setName("ScoreMatrixTest");
    f.setDescription("Compute an all-paths alignment between a read and template. For PPM visualization, try:\n"
    + "  " + f.getName() + " --ppm [...] | pnmscale 5 - | display -\n");
    f.registerOptional('t', TEMPLATE_FLAG, String.class, "DNA", "template bases", "AACCGGTTAACCGGTTTTTTGACTGACTGACGAGTACATGCTGCAGCGTAGCG");
    f.registerOptional('r', READ_FLAG, String.class, "DNA", "read bases",                "TAggCGGTTTTT GACTGACTcACGAtTtCATGCTGCAGaGT");
    f.registerOptional('q', QUALITY_FLAG, Integer.class, CommonFlags.INT, "default phread quality score", 20);
    f.registerOptional('s', START_FLAG, Integer.class, CommonFlags.INT, "start position of read", 7);
    f.registerOptional('m', MAXSHIFT_FLAG, Integer.class, CommonFlags.INT, "maxshift", 7);
    f.registerOptional('i', ITERATE_FLAG, Integer.class, CommonFlags.INT, "perform timing comparison using this many iterations");
    f.registerOptional('p', PPM_FLAG, "Output as PPM image");
    f.registerOptional('a', ASCII_FLAG, "Output as ascii matrix");
    if (!f.setFlags(args)) {
      System.exit(1);
    }
    final String tmpl = ((String) f.getValue(TEMPLATE_FLAG)).replaceAll(" ", "");
    final int start = (Integer) f.getValue(START_FLAG);
    final int maxShift = (Integer) f.getValue(MAXSHIFT_FLAG);
    final String read = ((String) f.getValue(READ_FLAG)).replaceAll(" ", "");
    final double[] quality = new double[read.length()];
    Arrays.fill(quality, MathUtils.phredToProb((Integer) f.getValue(QUALITY_FLAG)));
    final Environment env = new EnvironmentImplementation(
      maxShift,
      DNA.stringDNAtoByte(tmpl),
      start,
      DNA.stringDNAtoByte(read),
      quality
    );
    Exam.integrity(env);
    System.err.println("Template: " + tmpl);
    System.err.println("Read:     " + read);
    if (f.isSet(ITERATE_FLAG)) {
      final int repeat = (Integer) f.getValue(ITERATE_FLAG);
      for (PossibilityArithmetic arith : new PossibilityArithmetic[]{LogPossibility.SINGLETON, LogApproximatePossibility.SINGLETON, SimplePossibility.SINGLETON}) {
        final long time0 = System.nanoTime();
        ScoreMatrix sm1 = null;
        for (int i = 0; i < repeat; ++i) {
          sm1 = new ScoreMatrix(arith, new MockRealignParams());
          sm1.setEnv(env);
        }
        final long time1 = System.nanoTime();
        final double usec1 = (time1 - time0) / (1000.0 * repeat);
        System.out.println(StringUtils.padRight(arith.toString(), 30) + " time=" + usec1 + "usec  scoreLn=" + sm1.totalScoreLn());
      }
    } else {
      final AbstractAllPaths sm = new ScoreMatrix(LogApproximatePossibility.SINGLETON, RealignParamsGenome.SINGLETON);
      sm.setEnv(env);
      System.err.println("Probability: " + sm.totalScore());
      System.err.println("Log Score:   " + sm.totalScoreLn());
      if (f.isSet(PPM_FLAG)) {
        System.out.print(sm.toPpm());
      } else if (f.isSet(ASCII_FLAG)) {
        System.out.println(sm.toString());
      }
    }
  }
}
