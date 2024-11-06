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

package com.rtg.variant.util.arithmetic;

import junit.framework.TestCase;

/**
 * Common tests for Possibility implementations.
 */
public abstract class AbstractPossibilityTest extends TestCase {

  protected abstract PossibilityArithmetic arithmetic();

  protected abstract double tolerance();

  protected abstract String expectedToString();

  public void testGt() {
    final PossibilityArithmetic arith = arithmetic();
    assertTrue(arith.gt(arith.one(), arith.zero()));
    assertFalse(arith.gt(arith.zero(), arith.one()));
    checkGt(arith, 0.4, 0.5);
    checkGt(arith, 0.50, 0.5142);
    checkGt(arith, 0.51, 0.5142);
    checkGt(arith, 0.514, 0.5142);
    checkGt(arith, 0.5142, 0.5142);
    checkGt(arith, 0.1, 0.0);
    checkGt(arith, 0.01, 0.0);
    checkGt(arith, 0.001, 0.0);
    checkGt(arith, 0.0001, 0.0);
    checkGt(arith, 0.00001, 0.0);
    checkGt(arith, 0.000001, 0.0);
    checkGt(arith, 1.0, 0.9);
    checkGt(arith, 1.0, 0.99);
    checkGt(arith, 1.0, 0.999);
    checkGt(arith, 1.0, 0.9999);
    checkGt(arith, 1.0, 0.99999);
  }

  private void checkGt(final PossibilityArithmetic arith, final double a, final double b) {
    final double ap = arith.prob2Poss(a);
    final double bp = arith.prob2Poss(b);
    if (a > b + tolerance()) {
      assertTrue(arith.gt(ap, bp));
      assertFalse(arith.gt(bp, ap));
    } else if (b > a + tolerance()) {
      assertFalse(arith.gt(ap, bp));
      assertTrue(arith.gt(bp, ap));
    }
  }

  public void testIsZero() {
    final PossibilityArithmetic arith = arithmetic();
    assertTrue(arith.isZero(arith.zero()));
    assertFalse(arith.isZero(arith.one()));
  }

  public void testIsValid() {
    final PossibilityArithmetic arith = arithmetic();
    final double poss = arith.prob2Poss(0.0);
    assertTrue(arith.isValidPoss(poss));
    assertTrue(arith.isValidPoss(arith.prob2Poss(Double.MIN_VALUE)));
    assertTrue(arith.isValidPoss(arith.prob2Poss(0.1)));
    assertTrue(arith.isValidPoss(arith.prob2Poss(1.0)));

    assertFalse(arith.isValidPoss(arith.prob2Poss(1.1)));
    assertFalse(arith.isValidPoss(arith.prob2Poss(Double.MAX_VALUE)));
    assertFalse(arith.isValidPoss(arith.prob2Poss(Double.POSITIVE_INFINITY)));
    assertFalse(arith.isValidPoss(arith.prob2Poss(Double.NEGATIVE_INFINITY)));
    assertFalse(arith.isValidPoss(arith.prob2Poss(Double.NaN)));
  }

  public void testZero() {
    final PossibilityArithmetic arith = arithmetic();
    final double zero = arith.zero();
    assertEquals(0.0, arith.poss2Prob(zero));
    assertEquals(Double.NEGATIVE_INFINITY, arith.poss2Ln(zero));
    assertTrue(arith.isValidPoss(zero));
    assertEquals(zero, arith.prob2Poss(0.0));
    assertEquals(zero, arith.ln2Poss(Double.NEGATIVE_INFINITY));
  }

  public void testOne() {
    final PossibilityArithmetic arith = arithmetic();
    final double one = arith.one();
    assertEquals(1.0, arith.poss2Prob(one));
    assertEquals(0.0, arith.poss2Ln(one), 0.0);
    assertTrue(arith.isValidPoss(one));
    assertEquals(one, arith.prob2Poss(1.0));
    assertEquals(one, arith.ln2Poss(0.0));
  }

  private void checkAdd(final PossibilityArithmetic arith, final double p1, final double p2, double tolerance) {
    final double p = p1 + p2;
    final double s1 = arith.prob2Poss(p1);
    assertTrue(arith.isValidPoss(s1));
    final double s2 = arith.prob2Poss(p2);
    assertTrue("" + s2, arith.isValidPoss(s2));
    final double s = arith.add(s1, s2);
    assertTrue(arith.isValidPoss(s));
    final double pp = arith.poss2Prob(s);
    assertEquals(p, pp, tolerance);
  }

  public void testAdd() {
    final PossibilityArithmetic arith = arithmetic();
    checkAdd(arith, 0.0, 0.0, 0.0);
    checkAdd(arith, 0.0, 1.0, 0.0);
    checkAdd(arith, 1.0, 0.0, 0.0);
    checkAdd(arith, 0.25, 0.25, tolerance());
    checkAdd(arith, 0.1234, 0.678, tolerance());
    checkAdd(arith, Double.MIN_NORMAL, Double.MIN_NORMAL, tolerance());
    checkAdd(arith, Double.MIN_VALUE, Double.MIN_VALUE, tolerance());
  }

  private void checkMult(final PossibilityArithmetic arith, final double p1, final double p2, double tolerance) {
    final double p = p1 * p2;
    final double s1 = arith.prob2Poss(p1);
    assertTrue(arith.isValidPoss(s1));
    final double s2 = arith.prob2Poss(p2);
    assertTrue(arith.isValidPoss(s2));
    final double s = arith.multiply(s1, s2);
    assertTrue(arith.isValidPoss(s));
    final double pp = arith.poss2Prob(s);
    assertEquals(p, pp, tolerance);
  }

  public void testMult() {
    final PossibilityArithmetic arith = arithmetic();
    checkMult(arith, 0.0, 1.0, 0.0);
    checkMult(arith, 0.0, 0.0, 0.0);
    checkMult(arith, 1.0, 0.0, 0.0);
    checkMult(arith, 1.0, 1.0, 0.0);
    checkMult(arith, 0.1234, 0.678, tolerance());
    checkMult(arith, Double.MIN_NORMAL, Double.MIN_NORMAL, 0.0);
    checkMult(arith, Double.MIN_VALUE, Double.MIN_VALUE, 0.0);
    checkMult(arith, 1.0, Double.MIN_VALUE, tolerance());
  }

  private void checkDiv(final PossibilityArithmetic arith, final double p1, final double p2, double tolerance) {
    final double p = p1 / p2;
    final double s1 = arith.prob2Poss(p1);
    assertTrue(arith.isValidPoss(s1));
    final double s2 = arith.prob2Poss(p2);
    assertTrue("" + s2, arith.isValidPoss(s2));
    final double s = arith.divide(s1, s2);
    assertTrue(arith.isValidPoss(s));
    final double pp = arith.poss2Prob(s);
    assertEquals(p, pp, tolerance);
  }

  public void testDiv() {
    final PossibilityArithmetic arith = arithmetic();
    checkDiv(arith, 0.0, 1.0, 0.0);
    checkDiv(arith, 1.0, 1.0, 0.0);
    checkDiv(arith, 0.1234, 0.678, tolerance());
    checkDiv(arith, Double.MIN_NORMAL, Double.MIN_NORMAL, 0.0);
    checkDiv(arith, Double.MIN_VALUE, Double.MIN_VALUE, 0.0);
    checkDiv(arith, Double.MIN_NORMAL, 0.345, tolerance());
  }

  private void checkPow(final PossibilityArithmetic arith, final double p1, final double p2, double tolerance) {
    final double p = Math.pow(p1, p2);
    final double s1 = arith.prob2Poss(p1);
    assertTrue(arith.isValidPoss(s1));
    final double s = arith.pow(s1, p2);
    final double pp = arith.poss2Prob(s);
    assertEquals(p, pp, tolerance);
  }

  public void testPow() {
    final PossibilityArithmetic arith = arithmetic();
    checkPow(arith, 0.1, 0, 0);
    checkPow(arith, 1, 0, 0);
    checkPow(arith, 1, 1, 0);
    checkPow(arith, 0.5, 3, tolerance());
  }

  public void testUnderflow() {
    final PossibilityArithmetic arith = arithmetic();
    assertTrue(arith.underflow(arith.zero()));
    assertFalse(arith.underflow(arith.one()));
  }

  public void testToString() {
    assertEquals(expectedToString(), arithmetic().toString());
  }


  public void testConversion() {
    final PossibilityArithmetic arith = arithmetic();
    final double poss = arith.prob2Poss(0.5);
    assertEquals(poss, arith.poss2Poss(SimplePossibility.SINGLETON.prob2Poss(0.5), SimplePossibility.SINGLETON));
    assertEquals(poss, arith.poss2Poss(LogPossibility.SINGLETON.prob2Poss(0.5), LogPossibility.SINGLETON));
  }

}
