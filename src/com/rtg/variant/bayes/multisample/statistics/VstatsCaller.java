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

package com.rtg.variant.bayes.multisample.statistics;


import java.io.IOException;
import java.io.OutputStream;
import java.util.List;

import com.rtg.variant.Variant;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.AlleleStatistics;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.MultisampleJointCaller;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.StatisticsSnp;

/**
 */
public class VstatsCaller implements MultisampleJointCaller  {

  private final AlleleBalanceTable mDiploidTable = new AlleleBalanceTable();

  private final AlleleBalanceTable mHaploidTable = new AlleleBalanceTable();

  private final OutputStream mHaploidOutput;
  private final OutputStream mDiploidOutput;

  /**
   * @param params used to locate where the output is to go.
   * @throws IOException when opening output.
   */
  VstatsCaller(final VariantParams params) throws IOException {
    mDiploidOutput = params.outStream("diploid_statistics");
    mHaploidOutput = params.outStream("haploid_statistics");
  }

  @Override
  public <D extends Description, T extends HypothesesPrior<D>> Variant makeCall(String templateName, int position, int endPosition, byte[] ref, List<ModelInterface<?>> models, HaploidDiploidHypotheses<T> hypotheses) {
    assert models.size() == 1;

    final ModelInterface<?> model = models.get(0);
    if (model == null) {
      return null;
    }
    final int refNt = ref[position];
    if (refNt == 0) {
      return null;
    }
    final StatisticsSnp statistics = (StatisticsSnp) model.statistics();
    final AlleleStatistics<?> counts = statistics.counts();

    final int refCount = counts.count(refNt - 1);
    final int coverage = statistics.coverage();
    if (model.haploid()) {
      mHaploidTable.update(refCount, coverage);
    } else {
      mDiploidTable.update(refCount, coverage);
    }
    return null;
  }


  @Override
  public void endOfSequence() {
  }

  @Override
  public void close() throws IOException {
    try (OutputStream hapOut = mHaploidOutput; OutputStream dipOut = mDiploidOutput) {
      hapOut.write(mHaploidTable.toString().getBytes());
      dipOut.write(mDiploidTable.toString().getBytes());
    }
  }
}
