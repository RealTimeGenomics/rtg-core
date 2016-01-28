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

package com.rtg.variant.bayes;

import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.HypothesesCommon;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 * @author kurt
 */
public abstract class AbstractAlleleBalanceTest extends TestCase {
  int findIndexByName(String name) {
    final Description description = getDescriptionCommon();
    final Hypotheses<?> hypotheses = getDescriptionHypothesesCommon(description);
    for (int i = 0; i < hypotheses.size(); i++) {
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
        final double balance = balance(new int[]{i, j, 0, 0}, new double[]{0.0, 0.0, 0.0, 0.0}, findIndexByName("B:D"));
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
    private AlleleStatisticsInt mMockCounts;

    /**
     * @param description about which statistics are being collected.
     */
    public MockStatistics(Description description, int[] counts, double[] errors) {
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
    public MockAlleleStatisticsInt(int[] counts, double[] errors, Description description) {
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
