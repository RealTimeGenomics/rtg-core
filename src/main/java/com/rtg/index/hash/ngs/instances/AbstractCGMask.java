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
package com.rtg.index.hash.ngs.instances;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.ImplementHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.reader.CgUtils;

/**
 * Base class for masks operating on CG version 1 reads.
 */
public abstract class AbstractCGMask extends ImplementHashFunction {

  /**
   * An indicator for factories compatible with CG version 1 reads.
   */
  public interface CGHashFunctionFactory extends HashFunctionFactory { }

  /**
   * @param windowLength number of codes to be included in a window (used when hash called).
   * @param readCall used in subclasses to process results of read hits.
   * @param templateCall used in subclasses to process results of template hits.
   */
  public AbstractCGMask(final int windowLength, final ReadCall readCall, final TemplateCall templateCall) {
    super(CgUtils.CG_RAW_READ_LENGTH, windowLength, readCall, templateCall, LONG_BITS - (CgUtils.CG_RAW_READ_LENGTH + 7));
  }

  private static final long SCORE_MASK0 = 0b00000000000000000000000000000001111111111L; // (1L << 10) - 1;
  private static final long SCORE_MASK1 = 0b11111111111111111111111110000000000000000L; //((1L << 25) - 1) << 16;

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
    final long v1, v2;
    if (mTemplateCall.isReverse()) {
      v1 = mValuesR0;
      v2 = mValuesR1;
    } else {
      v1 = mValuesF0;
      v2 = mValuesF1;
    }
    return mMemScore.fastScore(bit0, bit1, cgAdjust(v1), cgAdjust(v2));
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
    return mMemScore.indelScore(bit0, bit1, cgAdjust(v1), cgAdjust(v2));
  }

}
