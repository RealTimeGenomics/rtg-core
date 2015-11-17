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

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.ImplementHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.reader.CgUtils;

public abstract class AbstractCG2Mask extends ImplementHashFunction {

  /**
   * An indicator for factories compatible with CG version 2 reads.
   */
  public interface CG2HashFunctionFactory extends HashFunctionFactory { }

  /**
   * @param windowLength number of codes to be included in a window (used when hash called).
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   */
  public AbstractCG2Mask(final int windowLength, final ReadCall readCall, final TemplateCall templateCall) {
    super(CgUtils.CG2_RAW_READ_LENGTH, windowLength, readCall, templateCall, LONG_BITS - (CgUtils.CG2_RAW_READ_LENGTH - 4));
  }

  private static final long SCORE_MASK0 = 0b00000001111111111111111111L;
  private static final long SCORE_MASK3 = 0b11111111110000000000000000L;
  private static final long SCORE_MASK4 = 0b01111111111000000000000000L;

  protected long cgAdjust3(final long x) {
    final long a = x & AbstractCG2Mask.SCORE_MASK0;
    final long b = (x & AbstractCG2Mask.SCORE_MASK3) << 3;
    return a | b;
  }

  protected long cgAdjust4(final long x) {
    final long a = x & AbstractCG2Mask.SCORE_MASK0;
    final long b = (x & AbstractCG2Mask.SCORE_MASK4) << 4;
    return a | b;
  }

  @Override
  public int fastScore(final int readId) {
    //System.err.println("fastScore readId=" + readId + " rc=" + mTemplateCall.isReverse());
    final long bit0 = mReadSequencesF1[readId];
    final long bit1 = mReadSequencesF2[readId];
    final long v1, v2;
    if (mTemplateCall.isReverse()) {
      v1 = mValuesR0;
      v2 = mValuesR1;
    } else {
      v1 = mValuesF0;
      v2 = mValuesF1;
    }
    final int f3 = mMemScore.fastScore(bit0, bit1, cgAdjust3(v1), cgAdjust3(v2));
    final int f4 = mMemScore.fastScore(bit0, bit1, cgAdjust4(v1), cgAdjust4(v2));
    return Math.min(f3, f4);
  }

  @Override
  public int indelScore(final int readId) {
    //System.err.println("indelScore rc=" + mTemplateCall.isReverse());
    final long bit0 = mReadSequencesF1[readId];
    final long bit1 = mReadSequencesF2[readId];
    final long v1, v2;
    if (mTemplateCall.isReverse()) {
      v1 = mValuesR0;
      v2 = mValuesR1;
    } else {
      v1 = mValuesF0;
      v2 = mValuesF1;
    }
    final int f3 = mMemScore.indelScore(bit0, bit1, cgAdjust3(v1), cgAdjust3(v2));
    final int f4 = mMemScore.indelScore(bit0, bit1, cgAdjust4(v1), cgAdjust4(v2));
    return Math.min(f3, f4);
  }
}
