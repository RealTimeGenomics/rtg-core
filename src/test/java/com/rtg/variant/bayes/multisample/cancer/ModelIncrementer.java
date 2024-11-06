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

package com.rtg.variant.bayes.multisample.cancer;

import java.util.List;

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Evidence;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.snp.EvidenceQ;

/**
 */
class ModelIncrementer<D extends Description> {
  private final List<ModelInterface<D>> mModels;
  private final double mDefaultQuality;
  private final double mMapQ;

  ModelIncrementer(List<ModelInterface<D>> models) {
    this(models, 0.05, 0.05);
  }

  ModelIncrementer(List<ModelInterface<D>> models, double defaultQuality) {
    this(models, defaultQuality, 0.05);
  }

  ModelIncrementer(List<ModelInterface<D>> models, double defaultQuality, double mapQ) {
    mModels = models;
    mDefaultQuality = defaultQuality;
    mMapQ = mapQ;
  }
  final ModelIncrementer<D> doRead(final int readNt, double quality, double mapQuality) {
    for (final ModelInterface<D> m : mModels) {
        final Evidence ev = new EvidenceQ(m.description(), readNt, 0, 0, mapQuality, quality, true, false, true, false, false);
        m.increment(ev);
    }
    return this;
  }

  final ModelIncrementer<D> doReads(final int readNt, double... qualities) {
    for (double quality : qualities) {
      doRead(readNt, quality, mMapQ);
    }
    return this;
  }

  final ModelIncrementer<D> doReads(final int numReads, final int readNt) {
    return doReads(numReads, readNt, mDefaultQuality);
  }

  final ModelIncrementer<D> doReadsCustom(final int numReads, Description desc, final int readNt, final double[] probabilities, final double pe, final double error, final double mapError) {
    for (int i = 0; i < numReads; ++i) {
      final Evidence ev = new EvidenceWrapper(desc, readNt, probabilities, error, mapError, pe);
      for (final ModelInterface<D> m : mModels) {
        m.increment(ev);
      }
    }
    return this;
  }

  final ModelIncrementer<D> doReads(final int numReads, final int readNt, final double quality) {
    for (int i = 0; i < numReads; ++i) {
      doRead(readNt, quality, mMapQ);
    }
    return this;
  }


  public List<ModelInterface<D>> freeze() {
    mModels.forEach(ModelInterface::freeze);
    return mModels;
  }
}
