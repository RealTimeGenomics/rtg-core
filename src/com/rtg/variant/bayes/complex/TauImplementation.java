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

import com.rtg.util.memo.DoubleFunction;
import com.rtg.util.memo.DoubleMemo;
import com.rtg.variant.dna.DNARangeNAT;
import com.rtg.variant.match.Match;

/**
 */
class TauImplementation {

  static double product(final Match match, final int matchStart, final String name, final int nameStart, final int nameEnd) {
    double prod = 1.0;
    for (int n = nameStart, m = matchStart; n < nameEnd; n++, m++) {
      final double q = match.baseError(m);
      final double s;
      if (q >= 0.75) { //allow for quality scores of 0 or 1 - prevent probability of 0
        s = 0.25;
      } else {
        final int nt = match.read(m);
        if (nt < 0) {
          s = 0.25;
        } else if (DNARangeNAT.DNA.toChar(nt) == name.charAt(n)) { //TODO check this is correct
          s = 1.0 - q;
        } else {
          s = q / 3.0;
        }
      }
      prod *= s;
      //System.err.println("product ms=" + matchStart + " name=" + name + " nameStart=" + nameStart + " nameEnd=" + nameEnd + " q=" + q + " s=" + s + " prod=" + prod);
    }
    return prod;
  }

  private final AbstractKappa mKappa;

  private final DoubleMemo mQMemo;

  /**
   * @param kappa used to build tau.
   */
  public TauImplementation(AbstractKappa kappa) {
    mKappa = kappa;
    mQMemo = new DoubleMemo(new DoubleFunction() {
      @Override
      public double fn(int i) {
        return qLocal(i);
      }
    });
  }

  /**
   * Probability that name can be transformed to match using probabilities from <code>kappaDistr</code>.
   * @param match evidence.
   * @param name hypothesis.
   * @return probability.
   */
  public double probability(final Match match, final String name) {
    if (match.isFixedLeft() && match.isFixedRight()) {
      return tau(match, name);
    } else if (!match.isFixedLeft() && match.isFixedRight()) {
      return taur(match, name);
    } else if (match.isFixedLeft() && !match.isFixedRight()) {
      return taul(match, name);
    } else if (!match.isFixedLeft() && !match.isFixedRight()) {
      return taub(match, name);
    } else {
      throw new RuntimeException();
    }
  }

  double taul(final Match match, final String name) {
    double sum = 0.0;
    final int l = name.length();
    for (int r = 0; r <= l; r++) {
      final String na = name.substring(0, r);
      final double tau = tau(match, na);
      sum += tau;
    }
    final double res = sum / (l + 1);
    //System.err.println("res=" + res);
    assert 0.0 <= res && res <= 1.0 && !Double.isNaN(res);
    return res;
  }

  double taur(final Match match, final String name) {
    double sum = 0.0;
    final int l = name.length();
    for (int r = 0; r <= l; r++) {
      final String na = name.substring(r);
      final double tau = tau(match, na);
      sum += tau;
      //System.err.println("m=" + match + " l=" + name + " na=" + na + " tau=" + tau + " wt=" + wt);
    }
    final double res = sum / (l + 1);
    //System.err.println("res=" + res);
    assert 0.0 <= res && res <= 1.0 && !Double.isNaN(res);
    return res;
  }

  double taub(final Match match, final String name) {
    double sum = 0.0;
    final int l = name.length();
    final double res;
    for (int r = 0; r <= l; r++) {
      for (int s = r; s <= l; s++) {
        final String na = name.substring(r, s);
        final double tau = tau(match, na);
        //System.err.println("m=" + match + " l=" + name + " na=" + na + " tau=" + tau + " wt=" + wt);
        sum += tau;
      }
    }
    res = 2 * sum / (l * (l + 1));
    //System.err.println("res=" + res);
    assert 0.0 <= res && res <= 1.0 && !Double.isNaN(res);
    return res;
  }

  /**
   * Probability that name can be transformed to match using probabilities from <code>kappaDistr</code>.
   * This assumes match is fixed at both ends.
   * @param match evidence.
   * @param name hypothesis.
   * @return probability.
   */
  double tau(final Match match, final String name) {
    final int l = name.length(); // hypothesis
    final int m = match.length(); // evidence
    final double kappa = mKappa.kappa(m, l);
    final double corr;
    double sum = 0.0;
    if (l >= m) {
      // hypothesis is longer or equal to evidence
      for (int d = 0; d <= m; d++) {
        final double p = product(match, 0, name, 0, d) * product(match, d, name, d + l - m, l);
        sum += p;
      }
      corr = 1.0 / (m + 1);
    } else {
      // hypothesis is shorter than evidence
      for (int d = 0; d <= l; d++) {
        final double p = product(match, 0, name, 0, d) * product(match, m - (l - d), name, d, l);
        sum += p;
      }
      corr = Math.pow(4.0, l - m) / (l + 1);
    }
    //System.err.println("sum= " + sum + " hypothesislength=" + l + " matchlength=" + m + " kappa=" + kappa + "corre=" + corr + " prob=" + prob);
    return sum * kappa * corr;
  }

  /**
   * Return the q correction factor (see theory).
   * @param l length of hypothesis.
   * @return the value of q.
   */
  public double q(final int l) {
    return mQMemo.fn(l);
  }

  /**
   * Return the q correction factor (see theory).
   * @param l length of hypothesis.
   * @return the value of q.
   */
  private double qLocal(final int l) {
    if (l < 0) {
      throw new IndexOutOfBoundsException(String.valueOf(l));
    }
    double q = 0;
    for (int m = 0; true; m++) {
      final double k = mKappa.kappa(m, l);
      final double delta = Math.pow(4.0, -m) * k * k;
      q += delta;
      if (m > l && delta < q * 0.00001) {
        break;
      }
    }
    return q;
  }
}
