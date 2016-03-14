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

package com.rtg.metagenomics;

import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractSolverTest extends TestCase {

  @Override
  protected void setUp() throws Exception {
    Diagnostic.setLogStream();
    super.setUp();
  }

  @Override
  protected void tearDown() throws Exception {
    super.tearDown();
  }

  abstract AbstractSolver getSolver(boolean verbose);

  AbstractSolver getSolver() {
    return getSolver(false);
  }

  public void testIsFinite() {
    assertFalse(AbstractSolver.isFinite(Double.POSITIVE_INFINITY));
    assertFalse(AbstractSolver.isFinite(Double.NEGATIVE_INFINITY));
    assertFalse(AbstractSolver.isFinite(Double.NaN));

    assertTrue(AbstractSolver.isFinite(Double.MAX_VALUE));
    assertTrue(AbstractSolver.isFinite(-Double.MAX_VALUE));
    assertTrue(AbstractSolver.isFinite(0.0));
    assertTrue(AbstractSolver.isFinite(-1.0));
    assertTrue(AbstractSolver.isFinite(1.0));
  }


  private static class StraightLine extends Line {
    private final double mConstant;
    private final double mMultiplier;
    StraightLine(final double constant, final double multiplier) {
      mConstant = constant;
      mMultiplier = multiplier;
    }

    @Override
    public double value(final double delta) {
      return mMultiplier * (mConstant + delta);
    }

    @Override
    public double[] values(double delta) {
      return new double[] {value(delta), mMultiplier};
    }

    @Override
    public int derivativeOrder() {
      return 1;
    }
  }

  public void testSimpleLine() {
    final AbstractSolver solver = getSolver();
    assertEquals(0.5, solver.solveLine(new StraightLine(-0.5, 1), 0.0, 1.0, 0.01), 1e-5);
    assertEquals(0.5, solver.solveLine(new StraightLine(-0.5, 1), 0.01), 1e-5);

    for (double k = 0; k < 10; k++) {
      assertEquals(k, solver.solveLine(new StraightLine(-k, 1.0), 1e-10), 1e-5);
    }
  }


  public void testLinearSmall() {
    final AbstractSolver solver = getSolver();
    for (double k = 0; k < 10; k++) {
      assertEquals(k, solver.solveLine(new StraightLine(-k, 0.5 / (k + 1)), 1e-10), 1e-5);
    }
  }

  public void testInfinitePos() {
    assertEquals(0.0, getSolver().solveLine(new StraightLine(Double.POSITIVE_INFINITY, 1.0), 1e-10), 1e-5);
  }

  public void testInfiniteNeg() {
    assertEquals(1.0, getSolver().solveLine(new StraightLine(Double.NEGATIVE_INFINITY, 1.0), 1e-10), 1e-5);
  }

  public void testNaN() {
    assertEquals(0.0, getSolver().solveLine(new StraightLine(Double.NaN, 1.0), 1e-10), 1e-5);
  }


  private static class ZeroLine extends Line {
    @Override
    public double value(final double delta) {
      if (delta < 2.7) {
        return -1.0;
      }
      return 0.0;
    }
  }

  public void testZeroLine() {
    final AbstractSolver solver = getSolver();
    final double res = solver.solveLine(new ZeroLine(), 0.01);
    assertTrue("res=" + res, res >= 2.7 && res <= 4.0);
  }

  private class ConstantLine extends Line {
    private final double mConstant;
    ConstantLine(final double constant) {
      mConstant = constant;
    }
    @Override
    public double value(final double delta) {
      return mConstant;
    }

    @Override
    public double[] values(double delta) {
      return new double[] {value(delta), 0};
    }

    @Override
    public int derivativeOrder() {
      return 1;
    }
  }

  public void testConstant() {
    assertEquals(0.0, getSolver().solveLine(new ConstantLine(0.0), 1e-10), 1e-5);
    assertEquals(1.0, getSolver().solveLine(new ConstantLine(-1.0), 1e-10), 1e-5);
    assertEquals(0.0, getSolver().solveLine(new ConstantLine(1.0), 1e-10), 1e-5);
  }

  public void testInfiniteCL() {
    assertEquals(0.0, getSolver().solveLine(new ConstantLine(Double.POSITIVE_INFINITY), 1e-10), 1e-5);
    assertEquals(1.0, getSolver().solveLine(new ConstantLine(Double.NEGATIVE_INFINITY), 1e-10), 1e-5);
  }

  public void testNaNCL() {
    assertEquals(0.0, getSolver().solveLine(new ConstantLine(Double.NaN), 1e-10), 1e-5);
  }

  private class AsymptoteLine extends Line {
    private final double mConstant;
    AsymptoteLine(final double constant) {
      mConstant = constant;
    }
    @Override
    public double value(final double delta) {
      return 1.5 - 1.0 / (delta + mConstant);
    }

    @Override
    public double[] values(double delta) {
      return new double[] {value(delta), 1.0 / ((delta + mConstant) * (delta + mConstant))};
    }

    @Override
    public int derivativeOrder() {
      return 1;
    }
  }

  public void testAsymptote() {
    //assertEquals(0.666667, getSolver().solveLine(new AsymptoteLine(0), 0.0, 1.0, 1e-10), 1e-5);
    assertEquals(0.666667, getSolver().solveLine(new AsymptoteLine(0), 1e-10), 1e-5);
    assertEquals(0.0, getSolver().solveLine(new AsymptoteLine(-1), 1e-10), 1e-5);
  }

  private class QuadraticLine extends Line {
    private final double mConstant;
    private final double mMultiplier;
    QuadraticLine(final double constant, final double multiplier) {
      mConstant = constant;
      mMultiplier = multiplier;
    }

    @Override
    public double value(double delta) {
      return mMultiplier * delta * delta + mConstant;
    }
    @Override
    public double[] values(double delta) {
      return new double[] {value(delta), 2.0 * mMultiplier * delta};
    }

    @Override
    public int derivativeOrder() {
      return 1;
    }
  }

  public void testQuadratic() {
    assertEquals(1.0, getSolver().solveLine(new QuadraticLine(-1.0, 1.0), 0.0, 2.0, 1e-10), 1e-5);
    assertEquals(1.0, getSolver().solveLine(new QuadraticLine(-1.0, 1.0), 1e-10), 1e-5);
    assertEquals(1.41421, getSolver().solveLine(new QuadraticLine(-1.0, 0.5), 1e-10), 1e-5);
  }

  private class PointsLine extends Line {
    private final double[] mPoints = {
        1.3354059041289173E-50, -246389.65940910578,
        1.471979106702489E-25, -246389.65940910578,
        3.315448323999644E-13, -246389.63463380933,
        2.0371414743438716E-7, -231168.551499933,
        0.0156251018571566, -2.4082890935593605E7,
        0.015635051329281908, -2.395512683161661E7,
        0.015750919427443177, -2.2310942500940815E7,
        0.01606031824089952, -1.647223456176196E7,
        0.016423057009132734, -6809394.008097261,
        0.01659278110624337, -1186384.7200078666,
        0.01662340410711451, -94204.89979365468,
        0.016625910442428293, -3750.331983089447,
        0.016626012165997343, -75.67932790517807,
        0.0166260142396049, -0.7697538435459137,
        0.016626014260803998, -0.003931701183319092,
        0.016626014260912553, -1.049041748046875E-5,
        0.016626014260912845, 0.0,
        0.016629614953103798, 130242.7581641674,
        0.016633215645295043, 260819.25658258796,
        0.016640417029786088, 522974.47822523117,
        0.016654819819967276, 1051301.6092852056,
        0.01668362747393721, 2124154.0588736236,
        0.016741344505446128, 4338786.705658495,
        0.016859284903777746, 9138404.628668368,
        0.017125788701312122, 2.138446329024011E7,
        0.017828520393491514, 6.3823127171896756E7,
        0.01959672254608351, 2.5501233256724402E8,
        0.023442525664723837, 1.5901995965831804E9,
        0.03125000000016577, 3.779628718686144E10,
        0.0625, 4.644727411832196E16,
        0.125, 2.092333191816366E29,
        0.25, 4.6126361027628015E54,
        0.5, 2.245117217184488E105
    };

    PointsLine() {

    }

    private double x(int i) {
      return mPoints[2 * i];
    }

    private double y(int i) {
      return mPoints[2 * i + 1];
    }

    private int size() {
      return mPoints.length / 2;
    }

    @Override
    public double value(double delta) {
      return values(delta)[0];
    }

    @Override
    public double[] values(double delta) {
      final double[] res = new double[2];
      //int blah = -1;

      if (delta <= x(0)) {
        res[0] =  y(0);
        res[1] = 0.0;
        //blah = -2;
      } else if (delta >= x(size() - 1)) {
        res[0] = y(size() - 1);
        res[1] = Double.MAX_VALUE;
        //blah = -3;
      } else {
        double px = x(0);
        double py = y(0);
        for (int i = 1; i < size(); i++) {
          final double x = x(i);
          final double y = y(i);
          if (delta <= x) {
            res[0] = py + (y - py) * (delta - px) / (x - px);
            res[1] = (y - py) / (x - px);
            //blah = i;
            break;
          }
          px = x;
          py = y;
        }
      }
      //System.err.println(blah + "\t" + delta + " --> " + res[0] + "\t" + res[1]);
      return res;
    }

    @Override
    public int derivativeOrder() {
      return 1;
    }
  }

  public void testPoints() {
    assertEquals(0.016626, getSolver().solveLine(new PointsLine(), 0.0, 1.0, 1e-10), 1e-5);
//    assertEquals(1.0, getSolver().solveLine(new QuadraticLine(-1.0, 1.0), 1e-10), 1e-5);
//    assertEquals(1.41421, getSolver().solveLine(new QuadraticLine(-1.0, 0.5), 1e-10), 1e-5);
  }

}
