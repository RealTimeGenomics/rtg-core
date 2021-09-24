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

package com.rtg.variant.util.arithmetic;

/**
 */
public class SimplePossibilityTest extends AbstractPossibilityTest {

  @Override
  protected PossibilityArithmetic arithmetic() {
    return SimplePossibility.SINGLETON;
  }

  @Override
  protected double tolerance() {
    return 1.0e-16;
  }

  public void testUnderflowTheHardWay() {
    final PossibilityArithmetic arith = arithmetic();
    double v = arith.one();
    final double m = arith.prob2Poss(0.0001);
    for (int i = 0; i < 77; ++i) {
      v = arith.multiply(v, m);
      if (arith.underflow(v)) {
        //System.err.println(i);
        return;
      }
    }
    fail();
  }

  @Override
  protected String expectedToString() {
    return "SimplePossibility";
  }
}
