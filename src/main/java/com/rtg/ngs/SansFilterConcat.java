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
