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
package com.rtg.variant.bayes.multisample.lineage;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import com.rtg.util.Pair;
import com.rtg.util.StringUtils;
import com.rtg.variant.bayes.AbstractGenotypeMeasure;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Implementation of a factor as a one-dimensional array.
 *
 */
public class DefaultFactor extends AbstractFactor implements ToDefaultFactor {

  /**
   * Return the given factor in the default implementation.  If the factor
   * is a default factor then it is returned unchanged.
   * @param f factor
   * @return default factor
   */
  public static DefaultFactor asDefault(final Factor f) {
    if (f instanceof ToDefaultFactor) {
      return ((ToDefaultFactor) f).asDefault();
    }
    final DefaultFactor def = new DefaultFactor(f.arithmetic(), f.scope());
    final Map<Variable, Integer> map = new HashMap<>(f.scope().size());
    for (int k = 0; k < def.mPoss.length; ++k) {
      def.mPoss[k] = f.p(def.getMap(map, k));
    }
    return def;
  }

  @Override
  public DefaultFactor asDefault() {
    return this;
  }

  /**
   * Return a normalized version of the factor
   * @param f original factor
   * @return factor where all entries sum to one
   */
  public static DefaultFactor asNormalized(final Factor f) {
    final DefaultFactor a = asDefault(f);
    double sum = a.arithmetic().zero();
    for (final double p : a.mPoss) {
      sum = a.arithmetic().add(sum, p);
    }
    if (a.arithmetic().isZero(sum)) {
      throw new IllegalArgumentException();
    }
    final double[] bPoss = new double[a.mPoss.length];
    for (int k = 0; k < a.mPoss.length; ++k) {
      bPoss[k] = a.arithmetic().divide(a.mPoss[k], sum);
    }
    return new DefaultFactor(a.arithmetic(), a.mScope, bPoss);
  }

  private final List<Variable> mScope;

  // Encode the multidimensional variables in a single dimensional array.  To do this
  // the array mStride[k] says how far to move in mPoss for an increment of variable k.
  private final int[] mStride;
  private final double[] mPoss; // stored in the arithmetic of this factor

  /**
   * Create a unit factor for the given scope.  That is, the factor that is everywhere one.
   * @param arith arithmetic
   * @param scope variables
   * @return the unit factor
   */
  static DefaultFactor unit(final PossibilityArithmetic arith, final Set<Variable> scope) {
    final DefaultFactor f = new DefaultFactor(arith, scope);
    Arrays.fill(f.mPoss, arith.one());
    return f;
  }

  /**
   * Factor.
   * @param arith arithmetic
   * @param scope variables
   */
  DefaultFactor(final PossibilityArithmetic arith, final Set<Variable> scope) {
    this(arith, new ArrayList<>(scope), (double[]) null);
  }

  /**
   * Factor.
   * @param arith arithmetic
   * @param scope variables
   * @param poss possibilities to assign in correct stride order and in correct arithmetic
   */
  DefaultFactor(final PossibilityArithmetic arith, final List<Variable> scope, final double... poss) {
    super(arith);
    mScope = scope;
    assert scope.size() == scope.size();
    int pos = 1;
    mStride = new int[scope.size()];
    for (int i = 0; i < scope.size(); ++i) {
      final Variable v = scope.get(i);
      mStride[i] = pos;
      pos *= v.size();
    }
    if (poss == null) {
      mPoss = new double[pos];
    } else {
      if (poss.length != pos) {
        throw new IllegalArgumentException("Length mismatch, should be " + pos + " but was " + poss.length);
      }
      mPoss = poss;
    }
  }

  private Map<Variable, Integer> getMap(final Map<Variable, Integer> map, final int key) {
    map.clear();
    for (int i = 0; i < mScope.size(); ++i) {
      final Variable variable = mScope.get(i);
      final int assignment = (key / mStride[i]) % variable.size();
      map.put(variable, assignment);
    }
    return map;
  }

  private int getKey(Map<Variable, Integer> values) {
    int key = 0;
    for (int i = 0; i < mScope.size(); ++i) {
      final Variable var = mScope.get(i);
      final Integer val = values.get(var);
      if (val == null) {
        throw new IllegalArgumentException("scope error");
      }
      if (val < 0 || val > var.size()) {
        throw new IllegalArgumentException(values.toString());
      }
      key += mStride[i] * val;
    }
    return key;
  }

  @Override
  public double p(Map<Variable, Integer> values) {
    checkScope(values);
    return mPoss[getKey(values)];
  }

  @Override
  public Set<Variable> scope() {
    return new HashSet<>(mScope);
  }

  private int stride(final Variable v) {
    final int index = mScope.indexOf(v);
    if (index >= 0) {
      return mStride[index];
    }
    return 0;
  }

  /**
   * Internal multiply.
   * @param phi the factor to multiply
   * @return a new factor which is the factor product of this and phi
   */
  private Factor multiplyInternal(DefaultFactor phi) {
    final PossibilityArithmetic arith = arithmetic();
    final Set<Variable> scope = new HashSet<>(scope());
    scope.addAll(phi.scope());
    final List<Variable> varList = new ArrayList<>(scope);
    final int[] assignment = new int[scope.size()];
    int j = 0;
    int k = 0;
    int factorSize = 1;
    for (final Variable v : varList) {
      factorSize *= v.size();
    }
    final double[] p = new double[factorSize];
    final int[] stride = new int[scope.size()];
    final int[] phiStride = new int[scope.size()];
    for (int l = 0; l < assignment.length; ++l) {
      stride[l] = stride(varList.get(l));
      phiStride[l] = phi.stride(varList.get(l));
    }
    for (int i = 0; i < p.length; ++i) {
      p[i] = arith.multiply(mPoss[j], arith.poss2Poss(phi.mPoss[k], phi.arithmetic()));
      for (int l = 0; l < assignment.length; ++l) {
        if (++assignment[l] == varList.get(l).size()) {
          j -= (varList.get(l).size() - 1) * stride[l];
          k -= (varList.get(l).size() - 1) * phiStride[l];
          assignment[l] = 0;
        } else {
          j += stride[l];
          k += phiStride[l];
          break;
        }
      }
    }
    return new DefaultFactor(arithmetic(), varList, p);
  }

  @Override
  public Factor marginal(final Set<Variable> variablesToKeep) {
    final PossibilityArithmetic arith = arithmetic();
    final List<Variable> varListRes = new ArrayList<>();
    final int[] resultStride = new int[mScope.size()];
    int factorSize = 1;
    for (int i = 0; i < mScope.size(); ++i) {
      final Variable v = mScope.get(i);
      if (variablesToKeep.contains(v)) {
        varListRes.add(v);
        resultStride[i] = factorSize;
        factorSize *= v.size();
      }
    }
    final double[] p = new double[factorSize];
    Arrays.fill(p, arith.zero());
    final int[] assignment = new int[mScope.size()];
    int j = 0;
    for (double mPos : mPoss) {
      p[j] = arith.add(p[j], mPos);
      for (int l = 0; l < assignment.length; ++l) {
        if (++assignment[l] == mScope.get(l).size()) {
          j -= (mScope.get(l).size() - 1) * resultStride[l];
          assignment[l] = 0;
        } else {
          j += resultStride[l];
          break;
        }
      }
    }
    return new DefaultFactor(arithmetic(), varListRes, p);
  }

  @Override
  public Factor condition(Map<Variable, Integer> conditions) {
    final List<Variable> varListRes = new ArrayList<>(mScope);
    varListRes.removeAll(conditions.keySet());
    if (varListRes.size() != mScope.size() - conditions.size()) {
      throw new IllegalArgumentException("scope error");
    }
    int factorSize = 1;
    for (final Variable v : varListRes) {
      factorSize *= v.size();
    }
    int j = 0;
    for (Map.Entry<Variable, Integer> entry : conditions.entrySet()) {
      j += entry.getValue() * stride(entry.getKey());
    }
    final double[] p = new double[factorSize];
    final int[] assignment = new int[varListRes.size()];
    for (int i = 0; i < p.length; ++i) {
      p[i] = mPoss[j];
      for (int l = 0; l < assignment.length; ++l) {
        if (++assignment[l] == varListRes.get(l).size()) {
          j -= (varListRes.get(l).size() - 1) * stride(varListRes.get(l));
          assignment[l] = 0;
        } else {
          j += stride(varListRes.get(l));
          break;
        }
      }
    }
    return new DefaultFactor(arithmetic(), varListRes, p);
  }

  @Override
  public Factor multiply(Factor other) {
    return multiplyInternal(DefaultFactor.asDefault(other));
  }

  /**
   *
   * @return The variable assignments that have the highest score and odds in the arithmetic space for this factor
   */
  public Pair<Map<Variable, Integer>, Double> best() {
    int best = 0;
    double bestScore = arithmetic().zero();
    double rest = arithmetic().zero();
    for (int i = 0; i < mPoss.length; ++i) {
      if (arithmetic().gt(mPoss[i], bestScore)) {
        best = i;
        rest = arithmetic().add(rest, bestScore);
        bestScore = mPoss[i];
      } else {
        rest = arithmetic().add(rest, mPoss[i]);
      }
    }
    return new Pair<>(getMap(new HashMap<Variable, Integer>(), best), arithmetic().divide(bestScore, rest));
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    final Collection<Variable> scope = new TreeSet<>(scope());
    sb.append(StringUtils.join("\t", scope)).append("\t").append("Value").append(StringUtils.LS);
    final Map<Variable, Integer> map = new HashMap<>(scope.size());
    for (int i = 0; i < mPoss.length; ++i) {
      getMap(map, i);
      for (final Variable v : scope) {
        sb.append(map.get(v)).append("\t");
      }
      sb.append(mPoss[i]).append(StringUtils.LS);
    }
    return sb.toString();
  }

  /**
   * A measure backed by a factor.
   * Only valid for factors with scope size == 1
   */
  public static class FactorGenotypeMeasure extends AbstractGenotypeMeasure {
    private final DefaultFactor mMeasure;

    /**
     * @param f the factor that backs this measure
     * @param hypotheses the hypothesis space this measure covers
     */
    public FactorGenotypeMeasure(DefaultFactor f, Hypotheses<?> hypotheses) {
      super(hypotheses);
      assert f.scope().size() == 1;
      mMeasure = f;
    }

    @Override
    public PossibilityArithmetic arithmetic() {
      return mMeasure.arithmetic();
    }

    @Override
    public double measure(int hypothesis) {
      return mMeasure.mPoss[hypothesis];
    }

    @Override
    public int size() {
      return mMeasure.mPoss.length;
    }
  }
}
