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
package com.rtg.ngs;

import java.io.IOException;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.pairedend.ReadStatusListener;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Concatenates and filters unmated SAM files, filtering on absolutely nothing.
 */
@TestClass("com.rtg.ngs.UnfilteredPairedEndOutputProcessorTest")
class SansFilterConcat extends AbstractMulticoreFilterConcat {
  private final long mReadIdOffset;
  private final ReadStatusListener mListener;

  SansFilterConcat(NgsParams params, ReadStatusListener listener, long readIdOffset) {
    super(params);
    mReadIdOffset = readIdOffset;
    mListener = listener;
  }

  @Override
  protected AbstractSamResultsFilter makeFilter() {
    final String readGroupId = mParams.outputParams().readGroup() != null ? mParams.outputParams().readGroup().getReadGroupId() : null;
    final SamResultsSansFilter filter = new SamResultsSansFilter(mListener, mReadIdOffset, mParams.buildFirstParams().reader().copy(), mParams.buildSecondParams() != null ? mParams.buildSecondParams().reader().copy() : null, readGroupId, mParams.legacyCigars());
    filter.setBamOutput(mParams.outputParams().bam());
    if (mParams.outputParams().outputReadNames()) {
      try {
        filter.setReadNames(mParams.buildFirstParams().reader().names());
      } catch (final IOException e) {
        filter.setReadNames(null);
        Diagnostic.warning("Failed to retrieve read names, using read id instead.");
      }
    }
    return filter;
  }
}
