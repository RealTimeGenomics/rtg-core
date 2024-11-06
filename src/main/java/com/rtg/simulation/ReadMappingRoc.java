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

package com.rtg.simulation;

import static com.rtg.vcf.eval.RocContainer.RocColumns.FALSE_POSITIVES;
import static com.rtg.vcf.eval.RocContainer.RocColumns.SCORE;
import static com.rtg.vcf.eval.RocContainer.RocColumns.TRUE_POSITIVES;

import java.text.DecimalFormat;
import java.util.Arrays;

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
  private final String mLabel;

  /**
   * @param label name of the attribute used for scoring
   */
  ReadMappingRoc(String label) {
    this(label, true);
  }

  /**
   * @param label name of the attribute used for scoring
   * @param ascending set this to true if low scores are good, otherwise set to false
   */
  ReadMappingRoc(String label, boolean ascending) {
    mAscending = ascending;
    mLabel = label;
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
    for (int i = 0; i <= getMaxScore(); ++i) {
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
    sb.append("#score field: ").append(mLabel).append(StringUtils.LS);
    sb.append("#");
    sb.append(String.join("\t", Arrays.asList(SCORE, TRUE_POSITIVES, FALSE_POSITIVES)));
    sb.append(StringUtils.LS);
    sb.append("255\t0.00\t0.00").append(StringUtils.LS);
    double tpCum = 0;
    double fpCum = 0;
    if (mAscending) {
      for (int i = 0; i <= getMaxScore(); ++i) {
        final double tpThis = getTp(i);
        final double fpThis = getFp(i);
        if (tpThis != 0 || fpThis != 0) {
          tpCum += tpThis;
          fpCum += fpThis;
          sb.append(i).append("\t").append(mFormatter.format(tpCum)).append("\t").append(mFormatter.format(fpCum)).append(StringUtils.LS);
        }
      }
    } else {
      for (int i = getMaxScore(); i >= 0; --i) {
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
