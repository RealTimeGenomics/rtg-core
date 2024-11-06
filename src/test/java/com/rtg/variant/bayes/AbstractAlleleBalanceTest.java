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

import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.HypothesesCommon;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractAlleleBalanceTest extends TestCase {

  int findIndexByName(String allele1, String allele2) {
    return findIndexByName(allele1 + VariantUtils.COLON + allele2);
  }
  int findIndexByName(String name) {
    final Description description = getDescriptionCommon();
    final Hypotheses<?> hypotheses = getDescriptionHypothesesCommon(description);
    for (int i = 0; i < hypotheses.size(); ++i) {
      if (hypotheses.name(i).equals(name)) {
        return i;
      }
    }
    throw new IllegalArgumentException("Unknown hypthesis name: " + name + " hypotheses: "  + hypotheses);
  }

  private DescriptionCommon getDescriptionCommon() {
    return new DescriptionCommon("B", "D", "E", "F");
  }

  public double balance(int[] counts, double[] errors, int hypothesis) {
    assert counts.length == errors.length;
    final AlleleBalanceProbability binomialAlleleBalance = getAlleleBalanceProbability();
    final double expected = 0.5;
    final Description description = getDescriptionCommon();
    final Statistics<?> statistics = new MockStatistics(description, counts, errors);
    final Hypotheses<?> hypotheses = getDescriptionHypothesesCommon(description);
    return binomialAlleleBalance.alleleBalanceLn(hypothesis, hypotheses, statistics);
  }
  public void dumpPlot() {
    for (int i = 0; i < 200; i += 5) {
      for (int j = 0; j < 200; j += 5) {
        final double balance = balance(new int[]{i, j, 0, 0}, new double[]{0.0, 0.0, 0.0, 0.0}, findIndexByName("B", "D"));
        System.err.println(String.format("%d\t%d\t%f", i, j, balance));
      }
    }
  }

  protected abstract AlleleBalanceProbability getAlleleBalanceProbability();

  private HypothesesCommon<Description> getDescriptionHypothesesCommon(Description description) {
    return new HypothesesCommon<>(description, SimplePossibility.SINGLETON, false, 0);
  }

  class MockStatistics extends StatisticsSnp {

    private final int mCoverage;
    private final double mError;
    private final AlleleStatisticsInt mMockCounts;

    /**
     * @param description about which statistics are being collected.
     */
    MockStatistics(Description description, int[] counts, double[] errors) {
      super(description);
      int coverage = 0;
      for (int count : counts) {
        coverage += count;
      }
      mCoverage = coverage;
      double totalError = 0.0;
      for (double error : errors) {
        totalError += error;
      }
      mError = totalError;
      mMockCounts = new MockAlleleStatisticsInt(counts, errors, mDescription);
    }

    @Override
    public AlleleStatisticsInt counts() {
      return mMockCounts;
    }

    @Override
    protected void incrementBest(EvidenceInterface evidence, int bestHyp) {
      //TODO
    }

    @Override
    public int coverage() {
      return mCoverage;
    }

    @Override
    public double totalError() {
      return mError;
    }
  }

  private class MockAlleleStatisticsInt extends AlleleStatisticsInt {
    private final double[] mErrors;
    private final int[] mCounts;

    /**
     * Create a new empty count object.
     *
     * @param description set of possible values
     */
    MockAlleleStatisticsInt(int[] counts, double[] errors, Description description) {
      super(description);
      this.mCounts = counts;
      this.mErrors = errors;
    }

    @Override
    public double count(int index) {
      return mCounts[index];
    }

    @Override
    public double error(int index) {
      return mErrors[index];
    }
  }
}
