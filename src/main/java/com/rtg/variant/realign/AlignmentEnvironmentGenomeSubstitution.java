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

package com.rtg.variant.realign;

import com.rtg.variant.bayes.complex.ComplexTemplate;

/**
 * Replace a specified subsequence in a reference sequence.
 */
public final class AlignmentEnvironmentGenomeSubstitution implements AlignmentEnvironment {

  private final int mStart;
  private final int mEnd;

  private final int mTemplateLength;

  private final ComplexTemplate mCot;

  private final byte[] mReplace;

  private final int mReplacePosition;

  private final int mDelta;


  /**
   * Create an AlignmentEnvironment corresponding to a section of template, with some bytes replaced.
   * @param start position on template (0 based).
   * @param end position on template (0 based exclusive).
   * @param cot template with information about region being replaced.
   * @param replace the bytes that replace the region selected in cot.
   */
  public AlignmentEnvironmentGenomeSubstitution(int start, int end, ComplexTemplate cot, byte[] replace) {
    mStart = start;
    mEnd = end;
    mTemplateLength = cot.templateBytes().length;
    mCot = cot;
    mReplace = replace;
    mReplacePosition = cot.getStart() + replace.length;
    mDelta = replace.length - (cot.getEnd() - cot.getStart());
  }

  @Override
  public int start() {
    return mStart;
  }

  @Override
  public double quality(int index) {
    return RealignParamsGenome.SINGLETON.misMatch();
  }

  @Override
  public byte base(int index) {
    final int rx = index + mStart;
    final int ix;
    if (rx < mCot.getStart()) {
      ix = rx;
    } else if (rx >= mReplacePosition) {
      ix = rx - mDelta;
    } else {
      return mReplace[rx - mCot.getStart()];
    }
    if (ix < 0 || ix >= mCot.templateBytes().length) {
      return 0;
    }
    return mCot.templateBytes()[ix];
  }

  @Override
  public int subsequenceLength() {
    return mEnd - mStart + mDelta;
  }

  @Override
  public int templateLength() {
    return mTemplateLength;
  }

  @Override
  public boolean isInverted() {
    return false;
  }

}
