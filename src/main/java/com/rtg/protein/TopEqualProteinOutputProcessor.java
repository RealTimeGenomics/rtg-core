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
package com.rtg.protein;

import java.io.IOException;

import com.rtg.alignment.ActionsHelper;
import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.ngs.MapStatistics;
import com.rtg.ngs.NgsParams;
//import java.util.ArrayList;

/**
 * Thread safe implementation of Top Equal output processor for protein world
 */
public class TopEqualProteinOutputProcessor extends ProteinOutputProcessor {

  private final int mTopN;
  private final TopEqualProteinImplementation mTopEqual;

  static NgsParams checkTopNValue(final NgsParams params) {
    final int topN = params.outputParams().filter().topN();
    if (topN < 1 || topN > 250) {
      throw new IllegalArgumentException("Invalid topN");
    }
    return params;
  }

  /**
   * Construct a new master {@link TopEqualProteinOutputProcessor}
   * @param params {@link NgsParams} object
   * @param statistics collector of statistics.
   * @throws IOException if error
   */
  public TopEqualProteinOutputProcessor(NgsParams params, MapStatistics statistics) throws IOException {
    super(checkTopNValue(params), statistics);
    mTopN = params.outputParams().filter().topN();
    assert params.buildFirstParams().numberSequences() <= Integer.MAX_VALUE;
    final int numberSequences = (int) params.buildFirstParams().numberSequences();
    mTopEqual = new TopEqualProteinImplementation(mTopN, numberSequences);
  }

  private TopEqualProteinOutputProcessor(NgsParams params, TopEqualProteinOutputProcessor master, TopEqualProteinImplementation topEqual, SharedStatusCollector collector, MapStatistics statistics) throws IOException {
    super(checkTopNValue(params), master, 0, true, collector, statistics);
    mTopN = params.outputParams().filter().topN();
    mTopEqual = topEqual;
    mChildren = null;
  }

  @Override
  public void finish() throws IOException {
    if (this == mMaster) {
    closeChildren();
    for (int l = 0; l < mTopEqual.numResults(); ++l) {
      final int count = mTopEqual.resultCount(l);
      if (count > mTopN) {
        //exceeds counts
        mSharedStatusCollector.setStatus(l, SharedStatusCollector.EXCEEDS_N_THRESHOLD);
      }
      if (count > 0 && count <= mTopN) {
        for (int i = 0; i < count; ++i) {
          super.writeResult(mTopEqual.result((long) l * mTopN + i));
        }
      }
    }
    super.writeUnmapped();
    mSharedStatusCollector.calculateStatistics();
    }
  }

  @Override
      protected boolean retainResult(final int[] res, final int readId, int templateId, int readAndFrame, int plen) throws IOException {
      return super.retainResult(res, readId, templateId, readAndFrame, plen)
      && mTopEqual.score(readId) >= ActionsHelper.alignmentScore(res)
      && !mTopEqual.contains(templateId, ActionsHelper.zeroBasedTemplateStart(res), readId, readAndFrame);
  }

  @Override
  void writeResult(ProteinAlignmentResult res) {
    mTopEqual.insertResult(res);
  }

  @Override
  public OutputProcessor threadClone(HashingRegion region) throws IOException {
    if (region != HashingRegion.NONE) {
      throw new UnsupportedOperationException();
    }
    final TopEqualProteinOutputProcessor p = new TopEqualProteinOutputProcessor(super.getParams(), this, mTopEqual, mSharedStatusCollector, mStatistics);
    mChildren.add(p);
    return p;
  }

  @Override
  public void threadFinish() throws IOException {
    try {
      finish();
    } finally {
      close();
    }
  }
}
