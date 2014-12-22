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

      for (int l = 0; l < mTopn.numResults(); l++) {
        final int count = mTopn.resultCount(l);
        final int score = mTopn.edgeScore(l);
        int edgeCount = mTopn.edgeScoreCount(l);
        boolean written = false;
        for (long i = 0, j = (long) l * mN; i < count; i++, j++) {
          final ProteinAlignmentResult res = mTopn.result(j);
          if (res.alignmentScore() < score) {
            super.writeResult(res);
            written = true;
          } else if (mN - i >= edgeCount) {
            edgeCount--;
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
