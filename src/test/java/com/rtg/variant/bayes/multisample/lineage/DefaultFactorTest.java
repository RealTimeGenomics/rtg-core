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
package com.rtg.variant.bayes.multisample.lineage;

import static com.rtg.util.StringUtils.LS;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.rtg.util.Pair;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 * Test class.
 */
public class DefaultFactorTest extends TestCase {

  static List<Variable> makeScope(final Variable... var) {
    return Arrays.asList(var);
  }

  public void testP() {
    final Variable a = new Variable("a", 3);
    final Variable b = new Variable("b", 2);
    final List<Variable> scopePhi1 = makeScope(a, b);
    final double[] p = {1, 2, 3, 4, 5, 6};
    final DefaultFactor phi = new DefaultFactor(SimplePossibility.SINGLETON, scopePhi1, p);
    final Map<Variable, Integer> map = new HashMap<>();
    for (int bb = 0, k = 0; bb < b.size(); ++bb) {
      map.put(b, bb);
      for (int aa = 0; aa < a.size(); ++aa, ++k) {
        map.put(a, aa);
        assertEquals(p[k], phi.p(map), 1e-10);
      }
    }
    assertTrue(phi == DefaultFactor.asDefault(phi));
    assertEquals("a\tb\tValue" + LS + "0\t0\t1.0" + LS + "1\t0\t2.0" + LS + "2\t0\t3.0" + LS + "0\t1\t4.0" + LS + "1\t1\t5.0" + LS + "2\t1\t6.0" + LS, phi.toString());
    final Pair<Map<Variable, Integer>, Double> best = phi.best();
    assertEquals(6.0 / 15.0, best.getB(), 1e-5);
    final Map<Variable, Integer> resultMap = new HashMap<>();
    resultMap.put(a, 2);
    resultMap.put(b, 1);
    assertEquals(resultMap, best.getA());
    final Factor conditioned = phi.condition(Collections.singletonMap(a, 1));
    assertEquals("b\tValue" + LS + "0\t2.0" + LS + "1\t5.0" + LS, conditioned.toString());
  }

  public void testMultiply() {
    final Variable a = new Variable("a", 3);
    final Variable b = new Variable("b", 2);
    final Variable c = new Variable("c", 2);
    final List<Variable> scopePhi1 = makeScope(a, b);
    final List<Variable> scopePhi2 = makeScope(b, c);
    final double[] p1 = {1, 2, 3, 4, 5, 6};
    final double[] p2 = {0.5, 0.25, 0.125, 0.125};
    final DefaultFactor phi1 = new DefaultFactor(SimplePossibility.SINGLETON, scopePhi1, p1);
    final DefaultFactor phi2 = new DefaultFactor(SimplePossibility.SINGLETON, scopePhi2, p2);
    final Factor psi = phi1.multiply(phi2);
    final Set<Variable> scopePsi = psi.scope();
    assertEquals(new HashSet<>(makeScope(a, b, c)), scopePsi);
    assertEquals(SimplePossibility.SINGLETON, psi.arithmetic());
    final double[] expected = {0.5, 1, 1.5, 1, 1.25, 1.5, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75};
    final Map<Variable, Integer> map = new HashMap<>();
    for (int cc = 0, k = 0; cc < c.size(); ++cc) {
      map.put(c, cc);
      for (int bb = 0; bb < b.size(); ++bb) {
        map.put(b, bb);
        for (int aa = 0; aa < a.size(); ++aa, ++k) {
          map.put(a, aa);
          assertEquals("a=" + aa + " b=" + bb + " c=" + cc, expected[k], psi.p(map), 1e-10);
        }
      }
    }
    final Factor conditioned = psi.condition(Collections.singletonMap(a, 1));
    assertEquals("b\tc\tValue"
                 + LS + "0\t0\t1.0"
                 + LS + "1\t0\t1.25"
                 + LS + "0\t1\t0.25"
                 + LS + "1\t1\t0.625"
                 + LS
        , conditioned.toString());

    map.clear();
    map.put(a, 1);
    map.put(c, 1);
    final Factor conditioned2 = psi.condition(map);
    assertEquals("b\tValue"
                 + LS + "0\t0.25"
                 + LS + "1\t0.625"
                 + LS
        , conditioned2.toString());
//    System.err.println(psi);
  }

  public void testMarginal() {
    final Variable a = new Variable("a", 3);
    final Variable b = new Variable("b", 2);
    final List<Variable> scopePhi1 = makeScope(a, b);
    final double[] p = {1, 2, 3, 4, 5, 6};
    final DefaultFactor phi = new DefaultFactor(SimplePossibility.SINGLETON, scopePhi1, p);
    final Factor psi = phi.marginal(Collections.singleton(a));
    final Map<Variable, Integer> map = new HashMap<>();
    final double[] expected = {5, 7, 9};
    for (int aa = 0, k = 0; aa < a.size(); ++aa, ++k) {
      map.put(a, aa);
      assertEquals(expected[k], psi.p(map), 1e-10);
    }
  }

  public void testMarginalKeepAll() {
    final Variable a = new Variable("a", 3);
    final Variable b = new Variable("b", 2);
    final List<Variable> scopePhi1 = makeScope(a, b);
    final double[] p = {1, 2, 3, 4, 5, 6};
    final DefaultFactor phi = new DefaultFactor(SimplePossibility.SINGLETON, scopePhi1, p);
    final Factor psi = phi.marginal(new HashSet<>(scopePhi1));
    final Map<Variable, Integer> map = new HashMap<>();
    for (int bb = 0, k = 0; bb < b.size(); ++bb) {
      map.put(b, bb);
      for (int aa = 0; aa < a.size(); ++aa, ++k) {
        map.put(a, aa);
        assertEquals(p[k], psi.p(map), 1e-10);
      }
    }
  }

  public void testMarginalEmpty() {
    final Variable a = new Variable("a", 3);
    final Variable b = new Variable("b", 2);
    final List<Variable> scopePhi1 = makeScope(a, b);
    final double[] p = {1, 2, 3, 4, 5, 6};
    final DefaultFactor phi = new DefaultFactor(SimplePossibility.SINGLETON, scopePhi1, p);
    final Factor psi = phi.marginal(Collections.emptySet());
    assertEquals(21, psi.p(new HashMap<>()), 1e-10);
  }

  public void testMarginalNoSuchVariable() {
    final Variable a = new Variable("a", 3);
    final Variable b = new Variable("b", 2);
    final List<Variable> scopePhi1 = makeScope(a, b);
    final double[] p = {1, 2, 3, 4, 5, 6};
    final DefaultFactor phi = new DefaultFactor(SimplePossibility.SINGLETON, scopePhi1, p);
    try {
      phi.marginal(Collections.singleton(new Variable("c", 2)));
    } catch (final IllegalArgumentException e) {
      assertEquals("c not in scope", e.getMessage());
    }
  }

  public void testUnitFactor() {
    final Variable a = new Variable("a", 3);
    final Variable b = new Variable("b", 2);
    final HashSet<Variable> scope = new HashSet<>(makeScope(a, b));
    for (PossibilityArithmetic arith : new PossibilityArithmetic[] {SimplePossibility.SINGLETON, LogPossibility.SINGLETON}) {
      final DefaultFactor phi = DefaultFactor.unit(arith, scope);
      final Map<Variable, Integer> map = new HashMap<>();
      for (int bb = 0; bb < b.size(); ++bb) {
        map.put(b, bb);
        for (int aa = 0; aa < a.size(); ++aa) {
          map.put(a, aa);
          assertEquals(arith.one(), phi.p(map), 1e-10);
        }
      }
    }
  }

  public void testMeasure() {
      final DefaultFactor factor = new DefaultFactor(SimplePossibility.SINGLETON, Collections.singletonList(new Variable("G1", 5)), 0.4, 0.2, 0.1, 0.1, 0.2);
      final DefaultFactor.FactorGenotypeMeasure measure = new DefaultFactor.FactorGenotypeMeasure(factor, new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), false, 2));
      assertEquals(SimplePossibility.SINGLETON, measure.arithmetic());
      assertEquals(0, measure.best());
      assertEquals(2, measure.reference());
      assertEquals(5, measure.size());
      assertEquals(0.4 / 0.6, measure.bestPosterior());
      assertEquals(0.9 / 0.1, measure.nonIdentityPosterior());
  }
}
