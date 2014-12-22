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

import java.util.Arrays;
import java.util.List;

import com.reeltwo.jumble.annotations.JumbleIgnore;
import com.rtg.metagenomics.matrix.Matrix;
import com.rtg.metagenomics.matrix.Vector;
import com.rtg.util.SortedMultiSet;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class Frag extends IntegralAbstract {

  private final int[] mGenomes;

  private final int[] mCounts;

  private final int mTotalCount;

  private final int mN;

  // Number of copies of this instance
  private int mMultiplicity = 1;

  private static int chain(final int[] a, final int g) {
    assert g >= 0;
    int t = g;
    do {
      assert t >= 0;
      final int u = a[t];
      assert u <= t : u + ":" + t;
      if (u == -1 || u == t) {
        return t;
      }
      t = u;
    } while (true);
  }

  private static void setChain(final int[] a, final int g, final int v) {
    assert g >= 0;
    int t = g;
    do {
      assert t >= 0;
      final int u = a[t];
      assert u <= t : u + ":" + t;
      a[t] = v;
      if (u == -1 || u == t) {
        return;
      }
      t = u;
    } while (true);
  }

  Frag(int[] genomes, int[] counts, int totalCount, int n, int multiplicity) {
    mGenomes = genomes;
    mCounts = counts;
    mTotalCount = totalCount;
    mN = n;
    mMultiplicity = multiplicity;
  }

  /**
   * @param gens list of genomes contained in this fragment.
   */
  public Frag(final List<Integer> gens) {
    final int gn = gens.size();
    int last = -1;
    int n = 0;
    for (final int g : gens) {
      if (last != g) {
        n++;
      }
      last = g;
    }
    mGenomes = new int[n];
    mCounts = new int[n];
    mN = n;
    assert n > 0;

    int c = 1;
    int l = gens.get(0);
    int j = 0;
    for (int i = 1; ; i++) {
      if (i == gn || l != gens.get(i)) {
        //another one
        assert c > 0;
        mGenomes[j] = l;
        mCounts[j] = c;
        j++;
        c = 0;
      }
      if (i == gn) {
        break;
      }
      l = gens.get(i);
      c++;
    }
    mTotalCount = gens.size();
    assert j == n; // : "j=" + j + " n=" + n;
  }

  void init(final Vector mv) {
    final double m = multiplicity();
    for (int i = 0; i < mN; i++) {
      final int g = mGenomes[i];
      final int c = mCounts[i];
      final double incr = c * m / (double) mTotalCount;
      //m is limited by number of fragments (i.e. number of objects in RAM)
      //c is limited by the number of reads
      //therefore incr should be positive and not too large
      assert incr >= 0.0 && !Double.isInfinite(incr) && !Double.isNaN(incr) : incr;
      //the total in mv should also be ok as total of increments is limited by number of reads
      mv.incr(g, incr);
    }
  }

  /**
   * Compute Hessian in frequency space.
   * @param r frequency estimates.
   * @param hessian to be updated.
   * @return <code>ln(rho)</code>. (used for computing L).
   */
  double incrementR(final Vector r, final Matrix hessian) {
    double rho = 0.0;
    final double m = multiplicity();
    for (int i = 0; i < mN; i++) {
      final int g = mGenomes[i];
      final int c = mCounts[i];
      rho += r.get(g) * c;
    }
    //System.err.println("rho: " + rho);
    assert hessian.isSymmetric();
    for (int i = 0; i < mGenomes.length; i++) {
      final int gi = mGenomes[i];
      final int ci = mCounts[i];
      for (int j = i; j < mGenomes.length; j++) {
        final int gj = mGenomes[j];
        final int cj = mCounts[j];
        final double v = (m / (rho * rho)) * ci * cj;
        //System.err.println("i=" + i + " gi=" + gi + " j=" + j + " gj=" + gj + " m=" + m + " rho=" + rho + " ci=" + ci + " cj=" + cj + " v=" + v);
        hessian.incr(gi, gj, v);
      }
    }
    final double logRho = -Math.log(rho);
    //System.err.println("rho=" + rho + " ln(rho)=" + logRho);
    return logRho * m;
  }


  /**
   * Compute Jacobian and Hessian in log space.
   * @param r frequency estimates (NOT in log space).
   * @param jacobian to be updated (in log space).
   * @param hessian to be updated (in log space).
   * @return <code>ln(rho)</code>. (used for computing L).
   */
  double increment(final Vector r, final Vector jacobian, final Matrix hessian) {
    double rho = 0.0;
    final double m = multiplicity();
    for (int i = 0; i < mN; i++) {
      final int g = mGenomes[i];
      final int c = mCounts[i];
      rho += r.get(g) * c;
    }

    for (int i = 0; i < mN; i++) {
      final int g = mGenomes[i];
      final int c = mCounts[i];
      final double incr = -r.get(g) * c * m / rho;
      jacobian.incr(g, incr);
      hessian.incr(g, g, incr);
    }

    assert hessian.isSymmetric();
    for (int i = 0; i < mGenomes.length; i++) {
      final int gi = mGenomes[i];
      final int ci = mCounts[i];
      final double ri = r.get(gi) * ci / rho;
      for (int j = i; j < mGenomes.length; j++) {
        final int gj = mGenomes[j];
        final int cj = mCounts[j];
        final double rj = r.get(gj) * cj / rho;
        final double rij = ri * rj * m;
        //System.err.println("i=" + i + " gi=" + gi + " j=" + j + " gj=" + gj + " rij=" + rij);
        hessian.incr(gi, gj, rij);
      }
    }
    final double logRho = -Math.log(rho);
    //System.err.println("rho=" + rho + " ln(rho)=" + logRho);
    return logRho * m;
  }

  /**
   * Compute Jacobian.
   * Used by <code>SpeciesDirect</code> and uses the direct <code>r_i</code> values (not the logs).
   * @param r vector of current co-ordinates.
   * @param jacobian first partial derivatives in frequency space.
   * @return <code>ln(rho)</code>. (used for computing L).
   */
  double increment(final Vector r, final Vector jacobian) {
    double rho = 0.0;
    final double m = multiplicity();
    for (int i = 0; i < mN; i++) {
      final int g = mGenomes[i];
      final int c = mCounts[i];
      rho += r.get(g) * c;
    }

    for (int i = 0; i < mN; i++) {
      final int g = mGenomes[i];
      final int c = mCounts[i];
      final double incr = -c / rho;
      jacobian.incr(g, incr * m);
    }
    final double logRho = -Math.log(rho);
    //System.err.println("rho=" + rho + " ln(rho)=" + logRho);
    return logRho * m;
  }

  /**
   * Compute L.
   * Used by <code>SpeciesDirect</code> and uses the direct <code>r_i</code> values (not the logs).
   * @param r vector of current co-ordinates.
   * @return <code>ln(rho)</code>. (used for computing L).
   */
  double l(final Vector r) {
    double rho = 0.0;
    for (int i = 0; i < mN; i++) {
      final int g = mGenomes[i];
      final int c = mCounts[i];
      rho += r.get(g) * c;
    }
    final double logRho = -Math.log(rho);
    //System.err.println("rho=" + rho + " ln(rho)=" + logRho);
    return logRho * (double) multiplicity();
  }

  @Override
  public void toString(final StringBuilder sb) {
    for (int i = 0; i < mN; i++) {
      final int g = mGenomes[i];
      final int c = mCounts[i];
      sb.append("  ").append(g).append(":").append(c);
    }
    sb.append(" {").append(multiplicity()).append('}');
  }


  /**
   * @param a array to be updated with descending chains of connected genomes.
   */
  public void stats(final int[] a) {
    for (int i = 0; i < mN; i++) {
      final int gi = mGenomes[i];
      for (int j = 0; j < mN; j++) {
        final int gj = mGenomes[j];
        final int mi = chain(a, gi);
        final int mj = chain(a, gj);
        final int m = Math.min(mi, mj);
        setChain(a, gi, m);
        setChain(a, gj, m);
      }
    }
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mN > 0);
    Exam.assertEquals(mN, mGenomes.length);
    Exam.assertEquals(mN, mCounts.length);
    Exam.assertTrue(mTotalCount >= mN);
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    int last = 0;
    for (int i = 0; i < mN; i++) {
      final int g = mGenomes[i];
      Exam.assertTrue(g >= last);
      last = g;
      final int c = mCounts[i];
      Exam.assertTrue(c > 0);
    }
    return true;
  }

  void checkRange(final int range) {
    for (final int mGenome : mGenomes) {
      Exam.assertTrue(mGenome >= 0 && mGenome < range);
    }
  }

  /**
   * Compute a weighted sum of a vector against the fragment.
   * DOES NOT take account of multiplicity.
   * @param v vector used to weight the sum.
   * @return the sum.
   */
  public double sum(final Vector v) {
    double sum = 0.0;
    for (int i = 0; i < mN; i++) {
      final int g = mGenomes[i];
      final int c = mCounts[i];
      sum += v.get(g) * c;
    }
    return sum;
  }

  @Override
  public boolean equals(final Object obj) {
    if (obj instanceof Frag) {
      final Frag that = (Frag) obj;
      return Arrays.equals(that.mGenomes, mGenomes) && Arrays.equals(that.mCounts, mCounts);
    }
    return false;
  }

  @JumbleIgnore
  boolean identical(Frag that) {
    return this.equals(that)
        && this.mTotalCount == that.mTotalCount
        && this.mN == that.mN
        && this.mMultiplicity == that.mMultiplicity;
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(Arrays.hashCode(mGenomes), Arrays.hashCode(mCounts));
  }

  int multiplicity() {
    return mMultiplicity;
  }

  /**
   * @param mult multiplicity.
   */
  public void setMultiplicity(final int mult) {
    mMultiplicity = mult;
  }

  Frag subFrag(int... c) {
    final int[] newGenomes = new int[mGenomes.length];
    for (int i = 0; i < mGenomes.length; i++) {
      newGenomes[i] = c[mGenomes[i]];
    }
    return new Frag(newGenomes, mCounts, mTotalCount, mN, mMultiplicity);
  }

  int subBlock(int... d) {
    assert sameBlock(d);
    return d[mGenomes[0]];
  }

  @JumbleIgnore
  private boolean sameBlock(int[] d) {
    final int b = d[mGenomes[0]];
    for (int mGenome : mGenomes) {
      Exam.assertEquals(b, d[mGenome]);
    }
    return true;
  }

  /**
   * Update counts of the size of each fragment.
   * @param counts to be updated.
   */
  public void count(SortedMultiSet<Integer>[] counts) {
    final int cnt = mN;
    for (int i = 0; i < mN; i++) {
      counts[mGenomes[i]].add(cnt, mMultiplicity);
    }
  }
}
