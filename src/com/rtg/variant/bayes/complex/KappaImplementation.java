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

import java.util.Map;
import java.util.TreeMap;

import com.rtg.util.Utils;
import com.rtg.util.format.FormatInteger;
import com.rtg.util.format.FormatReal;
import com.rtg.util.integrity.Exam;
import com.rtg.variant.AbstractMachineErrorParams;

/**
 */
public class KappaImplementation  extends AbstractKappa {

  private final int mLo;

  private final double[] mPi;

  private final double[] mPiSum;

  private final double mIndelDecay;

  private final double mPiInfinityToStart;

  private final double mLastPi;

  private final double mPiToInfinity;

  /**
   * Simple constructor directly from parameters.
   * @param params used for construction.
   * @param indelLengthDecay rate at which indel priors decrease beyond the end of the explicit distributions.
   */
  //Used only for testing
  KappaImplementation(final AbstractMachineErrorParams params, final double indelLengthDecay) {
    this(
        params.errorInsEventRate(),
        params.errorInsDistribution(),
        params.errorDelEventRate(),
        params.errorDelDistribution(),
        indelLengthDecay
    );
  }

  /**
   * @param insEventRate event rate for insertions.
   * @param insDistribution distribution of insertion lengths. (0 length assumed 0.0)
   * @param delEventRate event rate for deletions.
   * @param delDistribution distribution of deletion lengths. (0 length assumed 0.0).
   * @param indelDecay rate at which indel priors decrease beyond the end of the explicit distributions.
   */
  KappaImplementation(final double insEventRate, final double[] insDistribution, final double delEventRate, final double[] delDistribution, final double indelDecay) {
    mIndelDecay = indelDecay;
    final Map<Integer, Double> kappaMap = new TreeMap<>();
    assert insDistribution[0] == 0.0;
    assert delDistribution[0] == 0.0;
    assert Exam.assertDistribution(insDistribution);
    assert Exam.assertDistribution(delDistribution);
    //Compute simple single parameter kappa
    putKappa0(insEventRate, delEventRate, kappaMap);
    putKappaInsert(insEventRate, insDistribution, kappaMap);
    putKappaDelta(delEventRate, delDistribution, kappaMap);

    mPi = createPi(kappaMap);
    mLo = -delDistribution.length + 1;
    mPiSum = createCumuPi(mPi);
    assert integrity();
    mPiInfinityToStart = mPi[0] * mIndelDecay  / (1.0 - mIndelDecay);
    assert mPiInfinityToStart >= 0.0 && !Double.isNaN(mPiInfinityToStart) && !Double.isInfinite(mPiInfinityToStart);
    mLastPi = mPi[mPi.length - 1];
    assert mLastPi >= 0.0 && !Double.isNaN(mLastPi) && !Double.isInfinite(mLastPi);
    mPiToInfinity = mPi[mPi.length - 1] * mIndelDecay / (1.0 - mIndelDecay);
    assert mPiToInfinity >= 0.0 && !Double.isNaN(mPiToInfinity) && !Double.isInfinite(mPiToInfinity);
  }

  private void putKappa0(final double insEventRate, final double delEventRate, final Map<Integer, Double> kappaMap) {
    final double k0 = 1.0 - insEventRate - delEventRate; //(1.0 - insEventRate) * (1.0 - delEventRate);
    assert k0 > 0.0 && !Double.isNaN(k0) && k0 <= 1.0;
    kappaMap.put(0, k0);
  }

  private void putKappaInsert(final double insEventRate, final double[] insDistribution, final Map<Integer, Double> kappaMap) {
    assert Exam.assertDistribution(insDistribution);
    assert insEventRate > 0.0 && insEventRate < 1.0 && !Double.isNaN(insEventRate) : insEventRate;

    for (int i = 1; i < insDistribution.length; ++i) {
      final double d = insEventRate * insDistribution[i];
      assert d >= 0.0 && !Double.isNaN(d) && d <= 1.0;
      kappaMap.put(i, d);
    }
  }

  private void putKappaDelta(final double delEventRate, final double[] delDistribution, final Map<Integer, Double> kappaMap) {
    assert Exam.assertDistribution(delDistribution);
    assert delEventRate > 0.0 && delEventRate < 1.0 && !Double.isNaN(delEventRate) : delEventRate;
    for (int i = 1; i < delDistribution.length; ++i) {
      final double d = delEventRate * delDistribution[i];
      assert d >= 0.0 && !Double.isNaN(d) && d <= 1.0;
      kappaMap.put(-i, d);
    }
  }

  static double[] createCumuPi(final double[] kappa) {
    final double[] cumu = new double[kappa.length];
    double c = 0.0;
    for (int i = kappa.length - 1; i >= 0; --i) {
      c += kappa[i];
      cumu[i] = c;
    }
    return cumu;
  }

  private static double[] createPi(final Map<Integer, Double> kappaMap) {
    final double[] pi = new double[kappaMap.size()];
    int k = 0;
    for (final Map.Entry<Integer, Double> e : kappaMap.entrySet()) {
      final double d = e.getValue();
      pi[k] = d;
      ++k;
    }
    return pi;
  }

  @Override
  double pi(final int i) {
    final int index = i - mLo;
    final double pi;
    if (index < 0) {
      pi = mPi[0] * Math.pow(mIndelDecay , -index);
    } else  if (index < mPi.length) {
      pi = mPi[index];
    } else {
      pi = mLastPi * Math.pow(mIndelDecay , index - mPi.length + 1);
    }
    assert pi >= 0.0 && !Double.isNaN(pi) && !Double.isInfinite(pi);
    return pi;
  }

  @Override
  double piSum(final int i) {
    final int index = i - mLo;
    final double piSum;
    if (index < 0) {
      final double corr1 = mPiInfinityToStart * (1.0 - Math.pow(mIndelDecay, -index));
      assert corr1 >= 0.0 && !Double.isNaN(corr1) && !Double.isInfinite(corr1);
      final double corr2 = mPiSum[0];
      assert corr2 >= 0.0 && !Double.isNaN(corr2) && !Double.isInfinite(corr2);
      piSum = corr1 + corr2 + mPiToInfinity;
    } else  if (index < mPi.length) {
      final double corr2 = mPiSum[index];
      assert corr2 >= 0.0 && !Double.isNaN(corr2) && !Double.isInfinite(corr2);
      piSum = corr2 + mPiToInfinity;
    } else {
      piSum = mPiToInfinity * Math.pow(mIndelDecay, index - mPi.length);
    }
    assert piSum >= 0.0 && !Double.isNaN(piSum) && !Double.isInfinite(piSum);
    return piSum;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("Kappa ").append(LS);
    sb.append("pi          ").append(Utils.realFormat(mPi, 6)).append(LS);
    sb.append("piSum       ").append(Utils.realFormat(mPiSum, 6)).append(LS);
    sb.append("corrections ").append(Utils.realFormat(mPiInfinityToStart, 6)).append("  ").append(Utils.realFormat(mLastPi, 6)).append("  ").append(Utils.realFormat(mPiToInfinity, 6)).append(LS);
    sb.append(" m\\l");
    final FormatInteger fi = new FormatInteger(10);
    final int hi = mPi.length - mLo + 1;
    for (int l = 0; l <= hi; ++l) {
      sb.append(fi.format(l));
    }
    sb.append(LS);
    final FormatInteger fm = new FormatInteger(4);
    final FormatReal fk = new FormatReal(1, 6);
    for (int m = 0; m <= hi; ++m) {
      sb.append(fm.format(m));
      for (int l = 0; l <= hi; ++l) {
        sb.append("  ").append(fk.format(kappa(m, l)));
      }
      sb.append(LS);
    }
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    final int lhi = mPi.length + mLo + 3;
    for (int l = 0; l <= lhi; ++l) {
      double sumk = 0.0;

      final int hi = mPi.length + mLo + l;
      for (int m = 0; m < hi; ++m) {
        sumk += kappa(m, l);
      }
      final double t = kappa(hi, l) / (1.0 - mIndelDecay);
      sumk += t;
      Exam.assertEquals(1.0, sumk, 0.00001);
    }
    return true;
  }

  @Override
  public boolean integrity() {
    //System.err.println(toString());
    Exam.assertTrue(mIndelDecay >= 0.0 && mIndelDecay < 1.0 && !Double.isNaN(mIndelDecay));
    for (final double p : mPi) {
      Exam.assertTrue(p >= 0.0 && !Double.isInfinite(p) && !Double.isNaN(p));
    }
    for (final double p : mPiSum) {
      Exam.assertTrue(p >= 0.0 && !Double.isInfinite(p) && !Double.isNaN(p));
    }
    Exam.assertTrue(mPiInfinityToStart >= 0.0 && !Double.isInfinite(mPiInfinityToStart) && !Double.isNaN(mPiInfinityToStart));
    Exam.assertTrue(mLastPi >= 0.0 && !Double.isInfinite(mLastPi) && !Double.isNaN(mLastPi));
    Exam.assertTrue(mLo <= -1);
    return true;
  }
}
