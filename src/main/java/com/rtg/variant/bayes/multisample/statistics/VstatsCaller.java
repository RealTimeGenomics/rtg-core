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

package com.rtg.variant.bayes.multisample.statistics;


import java.io.IOException;
import java.io.OutputStream;
import java.util.List;

import com.rtg.util.MathUtils;
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

    final int refCount = (int) MathUtils.round(counts.count(refNt - 1));
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
