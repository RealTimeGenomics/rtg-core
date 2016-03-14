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

package com.rtg.simulation;

import java.text.DecimalFormat;

import com.rtg.util.Histogram;
import com.rtg.util.StringUtils;

/**
 * Hold a simple ROC (integer scores, automatically combines values at the same score)
 */
class ReadMappingRoc {

  // Buys 2 decimal point of accuracy in weights
  private static final int SCALING = 100;

  private final DecimalFormat mFormatter = new DecimalFormat("0.00");
  private final Histogram mTpHistogram = new Histogram();
  private final Histogram mFpHistogram = new Histogram();
  private final boolean mAscending;

  ReadMappingRoc() {
    this(true);
  }

  /**
   * @param ascending set this to true if low scores are good, otherwise set to false
   */
  ReadMappingRoc(boolean ascending) {
    mAscending = ascending;
  }

  private void add(Histogram hist, int score, double weight) {
    hist.increment(score, (long) (weight * SCALING));
  }

  private double get(Histogram hist, int score) {
    return (double) hist.getValueUnbounded(score) / SCALING;
  }

  public void addTp(int score, double weight) {
    add(mTpHistogram, score, weight);
  }
  public void addFp(int score, double weight) {
    add(mFpHistogram, score, weight);
  }

  public void addTp(int score) {
    addTp(score, 1);
  }

  public void addFp(int score) {
    addFp(score, 1);
  }

  public double getTp(int score) {
    return get(mTpHistogram, score);
  }

  public double getFp(int score) {
    return get(mFpHistogram, score);
  }

  public int getMaxScore() {
    final int length = Math.max(mTpHistogram.getLength(), mFpHistogram.getLength());
    return length == 0 ? 0 : length - 1;
  }

  public String getDistribution() {
    final StringBuilder sb = new StringBuilder();
    sb.append("Score\ttrue_positives\tfalse_positives\terror_rate").append(StringUtils.LS);
    for (int i = 0; i <= getMaxScore(); i++) {
      final double tpThis = getTp(i);
      final double fpThis = getFp(i);
      if (tpThis != 0 || fpThis != 0) {
        sb.append(i).append("\t").append(mFormatter.format(tpThis)).append("\t").append(mFormatter.format(fpThis)).append("\t").append(fpThis / (tpThis + fpThis)).append(StringUtils.LS);
      }
    }
    return sb.toString();
  }

  // Something that can be plotted with our ROC tools
  public String getRoc(int total) {
    final StringBuilder sb = new StringBuilder();
    sb.append("#total baseline variants: ").append(total).append(StringUtils.LS);
    sb.append("#score field: MAPQ").append(StringUtils.LS);
    sb.append("#score\ttrue_positives\tfalse_positives").append(StringUtils.LS);
    sb.append("255\t0.00\t0.00").append(StringUtils.LS);
    double tpCum = 0;
    double fpCum = 0;
    if (mAscending) {
      for (int i = 0; i <= getMaxScore(); i++) {
        final double tpThis = getTp(i);
        final double fpThis = getFp(i);
        if (tpThis != 0 || fpThis != 0) {
          tpCum += tpThis;
          fpCum += fpThis;
          sb.append(i).append("\t").append(mFormatter.format(tpCum)).append("\t").append(mFormatter.format(fpCum)).append(StringUtils.LS);
        }
      }
    } else {
      for (int i = getMaxScore(); i >= 0; i--) {
        final double tpThis = getTp(i);
        final double fpThis = getFp(i);
        if (tpThis != 0 || fpThis != 0) {
          tpCum += tpThis;
          fpCum += fpThis;
          sb.append(i).append("\t").append(mFormatter.format(tpCum)).append("\t").append(mFormatter.format(fpCum)).append(StringUtils.LS);
        }
      }
    }
    return sb.toString();
  }

}
