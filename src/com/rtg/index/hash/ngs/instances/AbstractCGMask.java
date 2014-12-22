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
package com.rtg.index.hash.ngs.instances;

import com.rtg.index.hash.ngs.ImplementHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;

/**
 */
public abstract class AbstractCGMask extends ImplementHashFunction {

  /**
   * @param readLength number of nucleotides in a complete read.
   * @param windowSize number of codes to be included in a window (used when hash called).
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   */
  public AbstractCGMask(final int readLength, final int windowSize, final ReadCall readCall, final TemplateCall templateCall) {
    super(readLength, windowSize, readCall, templateCall);
  }

  /**
   * @param readLength number of nucleotides in a complete read.
   * @param windowLength number of codes to be included in a window (used when hash called).
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   * @param readLengthReverse how much to shift right when in reverse complement. Needs to be specified explicitly for CG masks.
   */
  public AbstractCGMask(final int readLength, final int windowLength, final ReadCall readCall, final TemplateCall templateCall, final int readLengthReverse) {
    super(readLength, windowLength, readCall, templateCall, readLengthReverse);
  }

  private static final long SCORE_MASK0 = (1L << 10) - 1;
  private static final long SCORE_MASK1 = ((1L << 25) - 1) << 16;

  protected long cgAdjust(final long x) {
    final long a = x & AbstractCGMask.SCORE_MASK0;
    final long b = (x & AbstractCGMask.SCORE_MASK1) >>> 6;
    return a | b;
  }

  @Override
  public int fastScore(final int readId) {
    //System.err.println("fastScore readId=" + readId + " rc=" + mTemplateCall.isReverse());
    final long bit0 = mReadSequencesF1[readId];
    final long bit1 = mReadSequencesF2[readId];
    if (mTemplateCall.isReverse()) {
      return mMemScore.fastScore(bit0, bit1, cgAdjust(mValuesR0), cgAdjust(mValuesR1));
    } else {
      return mMemScore.fastScore(bit0, bit1, cgAdjust(mValuesF0), cgAdjust(mValuesF1));
    }
  }

  @Override
  public int indelScore(final int readId) {
    //System.err.println("indelScore rc=" + mTemplateCall.isReverse());
    final long bit0 = mReadSequencesF1[readId];
    final long bit1 = mReadSequencesF2[readId];
    if (mTemplateCall.isReverse()) {
      return mMemScore.indelScore(bit0, bit1, cgAdjust(mValuesR0), cgAdjust(mValuesR1));
    } else {
      return mMemScore.indelScore(bit0, bit1, cgAdjust(mValuesF0), cgAdjust(mValuesF1));
    }
  }

}
