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
        final Evidence ev = new EvidenceQ(m.description(), readNt, 0, 0, mapQuality, quality, true, false, false, false);
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

  final ModelIncrementer<D> doReads(final int numReads, final int readNt, final double quality) {
    for (int i = 0; i < numReads; i++) {
      doRead(readNt, quality, mMapQ);
    }
    return this;
  }


  public List<ModelInterface<D>> freeze() {
    mModels.forEach(ModelInterface::freeze);
    return mModels;
  }
}
