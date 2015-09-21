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
