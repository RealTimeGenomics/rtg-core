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
      throw new IllegalArgumentException();
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
