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
package com.rtg.variant.bayes.complex;

import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.MachineErrorParams;
import com.rtg.variant.match.AlignmentMatch;
import com.rtg.variant.match.Match;

import junit.framework.TestCase;

/**
 */
public class TauImplementationTest extends TestCase {
  //see tau.xls for details

  private TauImplementation getMockTau() {
    final AbstractKappa kappa = new AbstractKappa() {
        private final double[][] mKappa = {
          {0.98, 0.101, 0.0202},
          {0.11, 0.97, 0.22},
          {0.02, 0.1, 0.95, 0.3},
          {0.0, 0.0, 0.002, 0.93, 0.004, 0.0008},
          {0.0, 0.0, 0.02, 0.03, 0.9, 0.05},
        };
        @Override
        public double kappa(int m, int l) {
          //System.err.println("m=" + m + " l=" + l);
          final int length = mKappa.length - 1;
          if (m > length) {
            //just for q testing
            final double d = mKappa[length][l];
            return d * Math.pow(0.2, m - length);
          }
          return mKappa[m][l];
        }
        @Override
        double pi(int i) {
          throw new UnsupportedOperationException();
        }
        @Override
        double piSum(int i) {
          throw new UnsupportedOperationException();
        }
        @Override
        public boolean integrity() {
          // do nothing
          return true;
        }
        @Override
        public void toString(StringBuilder sb) {
          // do nothing
        }
      };
    return new TauImplementation(kappa);
  }

  public void testProduct() {
    final Match m0 = new AlignmentMatch(null, "A", null, 0, 0, 1, 20);
    assertEquals(0.25, TauImplementation.product(m0, 0, "A", 0, 1));
    final Match m1 = new AlignmentMatch(null, "A", null, 1, 0, 1, 20);
    assertEquals(0.25, TauImplementation.product(m1, 0, "A", 0, 1));
    final Match m2 = new AlignmentMatch(null, "A", null, 2, 0, 1, 20);
    assertTrue(0.25 < TauImplementation.product(m2, 0, "A", 0, 1));
  }

  private void check(final double expected, final TauImplementation tau, final Match m, final String l) {
    assertEquals(expected, tau.tau(m, l), 0.000001);
    assertEquals(expected, tau.probability(m, l), 0.000001);
  }

  //test values used below in the left/right/both tests
  public void testTaulrb() {
    final TauImplementation tau = getMockTau();
    final Match m = new AlignmentMatch(null, "AC", null, 10, 0, 2, 20);
    check(0.769500, tau, m, "AC");
    check(0.011667, tau, m, "A");
    check(0.011667, tau, m, "C");
    check(0.001250, tau, m, "");

    final Match m1 = new AlignmentMatch(null, "A", null, 10, 0, 1, 20);
    check(0.102667, tau, m1, "AC");
    check(0.873000, tau, m1, "A");
    check(0.027500, tau, m1, "");
  }

  public void testTaul() {
    //see diploidcomplexscorertest.xls for details and testTau3
    final TauImplementation tau = getMockTau();
    final Match m1 = new AlignmentMatch(null, null, "AC", null, 10, 0, 2, 20, true, false);
    assertEquals(0.260806, tau.taul(m1, "AC"), 0.000001);
    assertEquals(0.260806, tau.probability(m1, "AC"), 0.000001);

    final Match m2 = new AlignmentMatch(null, null, "A", null, 10, 0, 1, 20, true, false);
    assertEquals(0.334389, tau.taul(m2, "AC"), 0.000001);
    assertEquals(0.334389, tau.probability(m2, "AC"), 0.000001);
  }

  public void testTaur() {
    //see diploidcomplexscorertest.xls for details and testTau3, testTaul
    final TauImplementation tau = getMockTau();
    final Match m1 = new AlignmentMatch(null, null, "AC", null, 10, 0, 2, 20, false, true);
    assertEquals(0.260806, tau.taur(m1, "AC"), 0.000001);
    assertEquals(0.260806, tau.probability(m1, "AC"), 0.000001);

    final Match m2 = new AlignmentMatch(null, null, "A", null, 10, 0, 1, 20, false, true);
    assertEquals(0.054167, tau.taur(m2, "AC"), 0.000001);
    assertEquals(0.054167, tau.probability(m2, "AC"), 0.000001);
  }

  public void testTaub() {
    //see diploidcomplexscorertest.xls for details and testTau3, testTaul
    final TauImplementation tau = getMockTau();
    final Match m1 = new AlignmentMatch(null, null, "AC", null, 10, 0, 2, 20, false, false);
    assertEquals(0.265528, tau.taub(m1, "AC"), 0.000001);
    assertEquals(0.265528, tau.probability(m1, "AC"), 0.000001);
  }

  public void testTau1() {
    //see diploidcomplexscorertest.xls for details
    final TauImplementation tau = getMockTau();
    final Match m1 = new AlignmentMatch(null, "ACGT", null, 10, 0, 4, 20);
    check(0.590490, tau, m1, "ACGT");
    check(0.013374, tau, m1, "ACGTT");
    check(0.001419, tau, m1, "ACG");
  }

  public void testTau2() {
    //see diploidcomplexscorertest.xls for details
    final TauImplementation tau = getMockTau();
    final Match m1 = new AlignmentMatch(null, "AAAA", null, 10, 0, 4, 20);
    check(0.590490, tau, m1, "AAAA");
    check(0.032805, tau, m1, "AAAAA");
    check(0.005468, tau, m1, "AAA");
  }

  public void testTau3() {
    final TauImplementation tau = getMockTau();
    final Match m1 = match("ACG");
    check(0.677970, tau, m1, "ACG");
    check(0.000757, tau, m1, "ACGT");
    check(0.000140, tau, m1, "AC");
  }

  private static TauImplementation defaultTau() {
    final AbstractMachineErrorParams params = MachineErrorParams.builder().create();
    return new TauImplementation(new KappaImplementation(params, 0.0));
  }
  //The following tests use the default priors and check only the ordering of matches
  public void testOrder() {
    checkOrder("ACCG", "ACCG", "ACCGT");
    checkOrder("ACCG", "ACCG", "ACC");
    checkOrder("ACCG", "ACCG", "ACG");
    checkOrder("ACCG", "ACCG", "CCG");
    checkOrder("ACCG", "ACTG", "ACTT");

    checkOrder("ACCG", "ACC", "ACT");
    checkOrder("ACCG", "ACC", "ATC");
    checkOrder("ACCG", "ACC", "TCC");
    checkOrder("ACCG", "ACC", "AC");
    checkOrder("ACCG", "ACC", "CC");

    checkOrder("CCG", "C", "A");
    checkOrder("CCG", "G", "A");

    checkOrder("", "G", "GG");
    checkOrder("", "GG", "GGG");
  }

  public void checkOrder(final String l, final String m1, final String m2) {
    final TauImplementation tau = defaultTau();
    final double score1 = tau.probability(match(m1), l);
    final double score2 = tau.probability(match(m2), l);
    assertTrue(" l=" + l + "m1=" + m1 + " score1=" + score1 + "m2=" + m2 + " score2=" + score2, score1 > score2);
  }

  public void testQ() {
    final TauImplementation tau = getMockTau();
    assertEquals(0.96345, tau.q(0), 0.00001);
    assertEquals(0.24606, tau.q(1), 0.00001);
    assertEquals(0.06892, tau.q(2), 0.00001);

    try {
      tau.q(-1);
    } catch (final IndexOutOfBoundsException e) {
      assertEquals("-1", e.getMessage());
    }
    //System.err.println("KappaImplementation.q = " + ip.q(0));

  }

  private Match match(final String m) {
    return new AlignmentMatch(null, m, null, 10, 0, m.length(), 20);
  }
}

