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

package com.rtg.variant.bayes.snp;

import com.rtg.variant.VariantOutputOptions;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.StatisticsInt;

/**
 * Maintains various counts etc.
 */
public class StatisticsSnp extends StatisticsInt {

  /**
   * @param description about which statistics are being collected.
   */
  public StatisticsSnp(final Description description) {
    super(description);
  }

  @Override
  public void addCountsToSample(VariantSample sample, ModelInterface<?> model, VariantOutputOptions params) {
    //only doing support statistics for snp calls
    final StringBuilder sb = new StringBuilder();
    mCounts.output(sb, '\t');
    sample.setStatisticsString(sb.toString());
    super.addCountsToSample(sample, model, params);
  }
}
