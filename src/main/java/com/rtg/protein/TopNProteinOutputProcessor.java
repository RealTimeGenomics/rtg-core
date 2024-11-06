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

/**
 * Thread safe implementation of Top Equal output processor for protein world
 */
public class TopNProteinOutputProcessor extends ProteinOutputProcessor {

  private final int mN;
  private final TopNProteinImplementation mTopn;

  /**
   * Construct a new {@link TopNProteinOutputProcessor}
   * @param params {@link NgsParams} object
   * @param statistics collector of statistics.
   * @throws IOException if error
   */
  public TopNProteinOutputProcessor(NgsParams params, MapStatistics statistics) throws IOException {
    super(TopEqualProteinOutputProcessor.checkTopNValue(params), statistics);
    mN = params.outputParams().filter().topN();
    assert params.buildFirstParams().numberSequences() <= Integer.MAX_VALUE;
    final int numberSequences = (int) params.buildFirstParams().numberSequences();
    mTopn = new TopNProteinImplementation(mN, numberSequences);
  }

  // Only used by children
  private TopNProteinOutputProcessor(NgsParams params, TopNProteinOutputProcessor master, TopNProteinImplementation topn, SharedStatusCollector collector) throws IOException {
    super(TopEqualProteinOutputProcessor.checkTopNValue(params), master, 0, true, collector, null);
    mN = params.outputParams().filter().topN();
    mChildren = null;
    mTopn = topn;
  }

  @Override
  public void finish() throws IOException {
    if (this == mMaster) {
      closeChildren();

      for (int l = 0; l < mTopn.numResults(); ++l) {
        final int count = mTopn.resultCount(l);
        final int score = mTopn.edgeScore(l);
        int edgeCount = mTopn.edgeScoreCount(l);
        boolean written = false;
        for (long i = 0, j = (long) l * mN; i < count; ++i, ++j) {
          final ProteinAlignmentResult res = mTopn.result(j);
          if (res.alignmentScore() < score) {
            super.writeResult(res);
            written = true;
          } else if (mN - i >= edgeCount) {
            --edgeCount;
            written = true;
            super.writeResult(res);
          }
        }
        if (!written && count != 0) {
          // we did have some results and we did not write any
          mSharedStatusCollector.setStatus(l, SharedStatusCollector.EXCEEDS_N_THRESHOLD);
        }
      }
      super.writeUnmapped();
      mSharedStatusCollector.calculateStatistics();
    }
  }

  @Override
      protected boolean retainResult(final int[] res, final int readId, int templateId, int readAndFrame, int plen) throws IOException {
      return super.retainResult(res, readId, templateId, readAndFrame, plen)
      && mTopn.edgeScore(readId) >= ActionsHelper.alignmentScore(res)
      && !mTopn.contains(templateId, ActionsHelper.zeroBasedTemplateStart(res), readId, readAndFrame);
  }

  @Override
  void writeResult(ProteinAlignmentResult res) {
    mTopn.insertResult(res);
  }

  @Override
  public OutputProcessor threadClone(HashingRegion region) throws IOException {
    if (region != HashingRegion.NONE) {
      throw new UnsupportedOperationException();
    }
    final TopNProteinOutputProcessor p = new TopNProteinOutputProcessor(super.getParams(), this, mTopn, mSharedStatusCollector);
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
