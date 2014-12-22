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

package com.rtg.variant.bayes.multisample.forwardbackward;


import java.util.List;

import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Factor;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.UnitFactor;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.family.MendelianAlleleProbability;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Calculations for forward backward algorithm,
 * see <code>pedigree.pdf</code> for more details especially of the notation used.
 * Currently only deals with monogamous families (each parent can only mate with at most one other parent).
 */
//TODO allow non-monogamous families.
final class CommonFormulas {
  private CommonFormulas() { }

  /**
   * See <code>pedigree.pdf</code> for more details especially of the notation used.
   * Note it is imperative that the values for the father and mother be given in the correct order.
   * @param eu <code>E(u)</code> values for the FATHER
   * @param ev <code>E(v)</code> values for the MOTHER
   * @param aIndex index into c of the current child.
   * @param c <code>C(j,k)_b</code> intermediate values from all siblings includes contributions from Mendelian table, S and B for children (in possibility space).
   * @param m Mendelian table.
   * @return the new A value.
   */
  static Factor<?> forwardA(final Factor<?> eu, final Factor<?> ev, final Hypotheses<?> ha, final int aIndex, final double[][][] c, final MendelianAlleleProbability m) {
    final Hypotheses<?> uhyp = eu.hypotheses();
    final int usize = eu.size();
    final Hypotheses<?> vhyp = ev.hypotheses();
    final int vsize = ev.size();
    final MutableFactor<?> a = new MutableFactor<>(ha, ha.arithmetic(), ha.size() == 0 ? 1 : ha.size());
    final Code code = CommonFormulas.maxCode(uhyp, vhyp);
    final int size = a.size();
    final PossibilityArithmetic arith = eu.arithmetic();
    for (int h = 0; h < size; h++) {
      double v = arith.zero();
      for (int j = 0; j < usize; j++) {
        final double euj = eu.p(j);
        for (int k = 0; k < vsize; k++) {
          final double evk = ev.p(k);
          final double mjkh = arith.ln2Poss(m.probabilityLn(code, j, k, h));
          final double ee = arith.multiply(euj, evk);
          double p = arith.multiply(ee, mjkh);
          for (int b = 0; b < c.length; b++) {
            if (b != aIndex) {
              p = arith.multiply(p, c[b][j][k]);
              //System.err.println("h=" + h + " j=" + j + " k=" + k + " v=" + String.format("%1.4f", arith.poss2Prob(p)));
            }
          }
          v = arith.add(v, p);
        }

      }
      a.set(h, v);
    }
    return a;
  }

  /**
   * @param a <code>A_u</code> above value.
   * @param s <code>S_u</code> singleton value.
   * @return <code>E(u)</code> intermediate value.
   */
  static Factor<?> forwardEMonogamous(final Factor<?> a, final Factor<?> s) {
    return dot(a, s);
  }

  /**
   * @param familyId id of family current E is being calculated for
   * @param a <code>A_u</code> above value.
   * @param s <code>S_u</code> singleton value.
   * @param b Bs for individual
   * @return <code>E(u)</code> intermediate value.
   */
  static Factor<?> forwardEGeneral(final int familyId, final Factor<?> a, final Factor<?> s, BContainer b) {
    Factor<?> ret = dot(a, s);
    for (int i = 0; i < b.size(); i++) {
      if (i != familyId) {
        ret = dot(ret, b.getB(i));
      }
    }
    return ret;
  }

  static Factor<?> computeD(Factor<?> s, BContainer b) {
    Factor<?> d = s;
    for (int i = 0; i < b.size(); i++) {
      d = CommonFormulas.dot(d, b.getB(i));
    }
    return d;
  }

  static Factor<?> computeQ(Factor<?> a, Factor<?> s, BContainer b) {
    return dot(a, computeD(s, b));
  }

  /**
   * @param fatherSize size of vectors for father.
   * @param motherSize size of vectors for mother.
   * @param code codes to use
   * @param d <code>D(b)</code> combination of model with below values.
   * @param m Mendelian table.
   * @return <code>C(j, k)_x</code>.
   */
  static double[][] backwardC(int fatherSize, int motherSize, Code code, final Factor<?> d, final MendelianAlleleProbability m) {
    //System.err.println("father=" + father);
    //System.err.println("mother=" + mother);
    //System.err.println("M=" + m);
    //final Hypotheses<?> hyp = d.hypotheses();
    final int size = d.size();
    final PossibilityArithmetic arith = d.arithmetic();
    final double[][] c;
    c = new double[fatherSize][motherSize];
    for (int j = 0; j < fatherSize; j++) {
      for (int k = 0; k < motherSize; k++) {
        double v = arith.zero();
        for (int h = 0; h < size; h++) {
          final double t = arith.multiply(d.p(h), arith.ln2Poss(m.probabilityLn(code, j, k, h)));
          v = arith.add(v, t);
        }
        c[j][k] = v;
      }
    }
    return c;
  }

  /**
   * Assumes that the hypotheses are none, haploid or diploid across the same base haploid set.
   * @param a first set of hypotheses.
   * @param b second set of hypotheses.
   * @return the most inclusive code covering the hypotheses.
   */
  static Code maxCode(final Hypotheses<?> a, final Hypotheses<?> b) {
    final int ap = a.ploidy().count();
    final int bp = b.ploidy().count();
    if (ap <= bp) {
      return b.code();
    } else {
      return a.code();
    }
  }

  /**
   * @param eu <code>E(u, v)</code> intermediate value for first parent, only needed for size.
   * @param ev <code>E(v, u)</code> intermediate value for other parent.
   * @param c <code>C(j,k)_b</code> indexed in order <code>[b - over all children][j - one parent's hypotheses][k - other parent's hypotheses]</code>
   * @param father true iff <code>eu</code> is the father
   * @return <code>B(u,v)</code> backward value for this parent.
   */
  static Factor<?> backwardB(final Factor<?> eu, final Factor<?> ev, final double[][][] c, boolean father) {
    //final Hypotheses<?> hyp = ev.hypotheses();
    final int size = ev.size();
    final int usize = eu.size();
    final PossibilityArithmetic arith = ev.arithmetic();
    assert ev.arithmetic() == eu.arithmetic();
    final MutableFactor<?> b = new MutableFactor<>(eu.hypotheses(), arith, usize);
    for (int j = 0; j < usize; j++) {
      double v = arith.zero();
      for (int k = 0; k < size; k++) {
        double p = ev.p(k);
        for (double[][] aC : c) {
          if (father) {
            p = arith.multiply(p, aC[j][k]);
          } else {
            p = arith.multiply(p, aC[k][j]);
          }
        }
        v = arith.add(v, p);
      }
      b.set(j, v);
    }
    return b;
  }

  /**
   * Compute dot product of two hypotheses vectors.
   * @param a first vector.
   * @param b second vector.
   * @return vector of products of values in the two vectors a and b.
   */
  static Factor<?> dot(final Factor<?> a, final Factor<?> b) {
    final Hypotheses<?> hyp = a.hypotheses();
    assert a.size() == b.size() : "a.size=" + a.size() + " b.size=" + b.size() + a.getClass().getName() + b.getClass().getName();
    assert a.arithmetic() == b.arithmetic() : a.arithmetic().getClass().getName() + " != " + b.arithmetic().getClass().getName();
    final PossibilityArithmetic arith = a.arithmetic();
    final MutableFactor<?> c = new MutableFactor<>(hyp, arith, a.size());
    //assert hyp == b.hypotheses();
    //final PossibilityArithmetic arith = hyp.arithmetic();
    //final MutableFactor<?> c = new MutableFactor<?>(hyp, a.size());
    for (int h = 0; h < a.size(); h++) {
      c.set(h, arith.multiply(a.p(h), b.p(h)));
    }
    return c;
  }

  /**
   * @param factor factor from which probabilities are to be extracted.
   * @return vector of values extracted from priors.
   */
  static Factor<?> createMutableFactor(final Factor<?> factor) {
    if (factor.hypotheses() instanceof HypothesesNone) {
      return new UnitFactor<>(factor.hypotheses(), factor.arithmetic(), 1);
    }
    return new MutableFactor<>(factor);
  }

  /**
   * @return initial values for A when no A calculations have been done on samples
   */
  static <T extends HypothesesPrior<?>> Factor<?>[] initialA(List<ModelInterface<?>> models, HaploidDiploidHypotheses<T> hypotheses) {
    final Factor<?>[] as = new Factor<?>[models.size()];
    for (int i = 0; i < models.size(); i++) {
      as[i] = CommonFormulas.createMutableFactor(hypotheses.get(models.get(i)));
    }
    return as;
  }

  static Factor<?> getUnit(Factor<?> model) {
    if (model.hypotheses() instanceof HypothesesNone) {
      return new UnitFactor<>(model.hypotheses(), model.hypotheses().arithmetic(), 1);
    }
    return new UnitFactor<>(model.hypotheses(), model.hypotheses().arithmetic());
  }

  /**
   * @return initial values for B when no B calculations have been done on samples
   */
  static BContainer[] initialB(List<ModelInterface<?>> models, int[] bsize) {
    final BContainer[] bs = new BContainer[models.size()];
    for (int i = 0; i < models.size(); i++) {
      final Factor<?> model = models.get(i);
      final Factor<?>[] init = new Factor<?>[bsize[i]];
      for (int j = 0; j < init.length; j++) {
        init[j] = getUnit(model);
      }
      bs[i] = new BContainer(init);
    }
    return bs;
  }

  /**
   * Create an <code>HypothesisVector</code> whose values are set from the priors in the <code>GenomePriorParams</code>.
   * @param params from which priors are taken.
   * @param hyp hypotheses underlying the vector.
   * @param arith used for representing values in the vector.
   * @return the <code>HypothesisVector</code> containing the priors.
   */
  public static Factor<?> prior2Vector(final GenomePriorParams params, final Hypotheses<?> hyp, final PossibilityArithmetic arith) {
    final int size = hyp.size();
    if (size == 0) {
      return new UnitFactor<>(hyp, hyp.arithmetic(), 1);
    }
    final MutableFactor<?> mhv = new MutableFactor<>(hyp, arith, size);
    final int ref = hyp.reference();
    double total = 0.0;
    for (int i = 0; i < size; i++) {
      final String name = hyp.name(i);
      final double prob = params.getPriorDistr(name)[ref];
      total += prob;
    }
    for (int i = 0; i < size; i++) {
      final String name = hyp.name(i);
      final double prob = params.getPriorDistr(name)[ref] / total;
      final double prior = arith.prob2Poss(prob);
      mhv.set(i, prior);
    }
    return mhv;
  }
}
