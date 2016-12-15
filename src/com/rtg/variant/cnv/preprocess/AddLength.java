/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.variant.cnv.preprocess;

import com.rtg.util.intervals.SequenceNameLocus;

/**
 * Computes region length
 */
class AddLength implements DatasetProcessor {
  @Override
  public void process(RegionDataset dataset) {
    final NumericColumn c = dataset.addColumn(new IntColumn("length"));
    for (SequenceNameLocus r : dataset.regions().values()) {
      c.add(r.getLength());
    }
  }
}
