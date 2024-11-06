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

package com.rtg.variant.bayes;

import static com.rtg.util.StringUtils.TAB;

import java.util.Random;

import com.rtg.util.TestUtils;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.util.VariantUtils;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractAlleleStatisticsTest<T extends AlleleStatistics<T>> extends TestCase {

  protected abstract T getAlleleStatistics(Description d);

  protected abstract void incrementAlleleStatistics(T stats, EvidenceInterface dist);


  public void test() {
    final T cn = getAlleleStatistics(DescriptionSnp.SINGLETON);
    for (int i = 0; i < 4; ++i) {
      assertEquals(0.0, cn.count(i));
      assertEquals(0.0, cn.error(i));
    }
    //System.err.println(cn.toString());
    TestUtils.containsAll(cn.toString(), " [0]  0  0.000 [1]  0  0.000 [2]  0  0.000 [3]  0  0.000");
  }

  public void test1() {
    final T cn = getAlleleStatistics(DescriptionSnp.SINGLETON);
    //System.err.println(cn.toString());
    final EvidenceInterface di = new EvidenceQ(DescriptionSnp.SINGLETON, 1, 0, 0, 0.0, VariantUtils.phredToProb(3), true, false, true, false, false);
    incrementAlleleStatistics(cn, di);
    TestUtils.containsAll(cn.toString(), " [0]  0  0.000 [1]  1  0.501 [2]  0  0.000 [3]  0  0.000");
  }

  EvidenceInterface di(final int read, final int score, double r) {
    return new EvidenceQ(DescriptionSnp.SINGLETON, read, 0, 0, r, VariantUtils.phredToProb(score), true, false, true, false, false);
  }

  public void test2() {
    final T cn = getAlleleStatistics(DescriptionSnp.SINGLETON);
    //System.err.println(cn.toString());
    incrementAlleleStatistics(cn, di(0, 3, 0.0));
    incrementAlleleStatistics(cn, di(1, 4, 0.0));
    incrementAlleleStatistics(cn, di(1, 4, 0.0));
    incrementAlleleStatistics(cn, di(2, 5, 0.0));
    incrementAlleleStatistics(cn, di(2, 5, 0.0));
    incrementAlleleStatistics(cn, di(2, 5, 0.0));
    incrementAlleleStatistics(cn, di(3, 6, 0.0));
    incrementAlleleStatistics(cn, di(3, 6, 0.0));
    incrementAlleleStatistics(cn, di(3, 6, 0.0));
    incrementAlleleStatistics(cn, di(3, 6, 0.0));
    TestUtils.containsAll(cn.toString(), " [0]  1  0.501 [1]  2  0.796 [2]  3  0.949 [3]  4  1.005");
    final StringBuilder sb = new StringBuilder();
    cn.output(sb, '\t');
    assertEquals(""
                         + TAB + "A" + TAB + "1" + TAB + "0.501"
                         + TAB + "C" + TAB + "2" + TAB + "0.796"
                         + TAB + "G" + TAB + "3" + TAB + "0.949"
                         + TAB + "T" + TAB + "4" + TAB + "1.005"
                        , sb.toString());
  }

  public void test3() {
    final T cn = getAlleleStatistics(DescriptionSnp.SINGLETON);
    //System.err.println(cn.toString());
    incrementAlleleStatistics(cn, di(0, 3, 0.1));
    incrementAlleleStatistics(cn, di(1, 4, 0.1));
    incrementAlleleStatistics(cn, di(1, 4, 0.1));
    incrementAlleleStatistics(cn, di(2, 5, 0.1));
    incrementAlleleStatistics(cn, di(2, 5, 0.1));
    incrementAlleleStatistics(cn, di(2, 5, 0.1));
    incrementAlleleStatistics(cn, di(3, 6, 0.1));
    incrementAlleleStatistics(cn, di(3, 6, 0.1));
    incrementAlleleStatistics(cn, di(3, 6, 0.1));
    incrementAlleleStatistics(cn, di(3, 6, 0.1));
    //see spreadsheet
//    assertEquals(-12.81030782 - EntropyRandom.RANDOM_PRIOR, cn.noisePosterior(), 1e-8);
    final String str = cn.toString();
    //System.err.println(str);
    TestUtils.containsAll(str, " [0]  1  0.551 [1]  2  0.917 [2]  3  1.154 [3]  4  1.304");
    final StringBuilder sb = new StringBuilder();
    cn.output(sb, '\t');
    assertEquals(""
                         + TAB + "A" + TAB + "1" + TAB + "0.551"
                         + TAB + "C" + TAB + "2" + TAB + "0.917"
                         + TAB + "G" + TAB + "3" + TAB + "1.154"
                         + TAB + "T" + TAB + "4" + TAB + "1.304"
                        , sb.toString());
  }

  public void testPerBaseError() {
    final T cn = getAlleleStatistics(DescriptionSnp.SINGLETON);
    //System.err.println(cn.toString());
    incrementAlleleStatistics(cn, di(0, 2, 0.0)); //this should be ignored in perBaseErrorCalculation
    incrementAlleleStatistics(cn, di(0, 3, 0.0));
    incrementAlleleStatistics(cn, di(1, 4, 0.0));
    incrementAlleleStatistics(cn, di(1, 4, 0.0));
    incrementAlleleStatistics(cn, di(2, 5, 0.0));
    incrementAlleleStatistics(cn, di(2, 5, 0.0));
    incrementAlleleStatistics(cn, di(2, 5, 0.0));
    incrementAlleleStatistics(cn, di(3, 6, 0.0));
    incrementAlleleStatistics(cn, di(3, 6, 0.0));
    incrementAlleleStatistics(cn, di(3, 6, 0.0));
    incrementAlleleStatistics(cn, di(3, 6, 0.0));
//    assertEquals(11, cn.coverage());


  }
  EvidenceInterface diPaired(final int read, final int score, double r, boolean mated) {
    return new EvidenceQ(DescriptionSnp.SINGLETON, read, 0, 0, r, VariantUtils.phredToProb(score), true, true, true, mated, false);
  }
  EvidenceInterface diStrand(final int read, final int score, double r, boolean forward) {
    return new EvidenceQ(DescriptionSnp.SINGLETON, read, 0, 0, r, VariantUtils.phredToProb(score), forward, true, true, true, false);
  }
  public void testAlleleBalance() {
    final T cn = getAlleleStatistics(DescriptionSnp.SINGLETON);
    final Random r = new Random();
    //System.err.println(cn.toString());
    for (int i = 0; i < 10; ++i) {
      incrementAlleleStatistics(cn, diStrand(0, 2, 0.0, r.nextBoolean()));
    }
    for (int i = 0; i < 6; ++i) {
      incrementAlleleStatistics(cn, diStrand(1, 2, 0.0, r.nextBoolean()));
    }
    assertEquals(2.1715, cn.alleleBalance(0, 1), 0.0001);
    assertEquals(19.5433, cn.alleleBalanceHomozygous(0, 16), 0.0001);
    incrementAlleleStatistics(cn, diStrand(3, 2, 0.0, r.nextBoolean()));
    assertEquals(25.0358, cn.alleleBalanceHomozygous(0, 17), 0.0001);
  }

  public void testUnmatedBias() {
    final T cn = getAlleleStatistics(DescriptionSnp.SINGLETON);
    //System.err.println(cn.toString());
    int mated = 0;
    int unmated = 0;
    for (int i = 0; i < 10; ++i) {
      incrementAlleleStatistics(cn, diPaired(0, 2, 0.0, true));
      ++mated;
    }
    for (int i = 0; i < 5; ++i) {
      incrementAlleleStatistics(cn, diPaired(0, 2, 0.0, false));
      ++unmated;
    }
    for (int i = 0; i < 4; ++i) {
      incrementAlleleStatistics(cn, diPaired(1, 2, 0.0, true));
      ++mated;
    }
    for (int i = 0; i < 11; ++i) {
      incrementAlleleStatistics(cn, diPaired(1, 2, 0.0, false));
      ++unmated;
    }
    final double prob = (double) unmated / (mated + unmated);
    assertEquals(5.2115, cn.unmatedBias(0, prob), 0.0001);
    assertEquals(5.2115, cn.unmatedBias(1, prob), 0.0001);
  }
  public void testStrandBias() {
    final T cn = getAlleleStatistics(DescriptionSnp.SINGLETON);
    //System.err.println(cn.toString());
    for (int i = 0; i < 10; ++i) {
      incrementAlleleStatistics(cn, diStrand(0, 2, 0.0, true));
    }
    for (int i = 0; i < 5; ++i) {
      incrementAlleleStatistics(cn, diStrand(0, 2, 0.0, false));
    }
    for (int i = 0; i < 4; ++i) {
      incrementAlleleStatistics(cn, diStrand(1, 2, 0.0, true));
    }
    for (int i = 0; i < 11; ++i) {
      incrementAlleleStatistics(cn, diStrand(1, 2, 0.0, false));
    }
    assertEquals(3.6191, cn.strandBias(0), 0.0001);
    assertEquals(7.0934, cn.strandBias(1), 0.0001);
  }
}
