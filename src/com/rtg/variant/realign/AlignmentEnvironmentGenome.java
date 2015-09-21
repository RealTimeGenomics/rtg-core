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


/**
 * Directly corresponds to an unmodified subsequence in a reference sequence.
 */
public final class AlignmentEnvironmentGenome implements AlignmentEnvironment {

  private final int mStart;
  private final int mEnd;

  private final byte[] mTemplate;

  /**
   * Create an AlignmentEnvironment corresponding to a section of template.
   * @param start position on template (0 based).
   * @param end position on template (0 based exclusive).
   * @param template supplies the full template bytes
   */
  public AlignmentEnvironmentGenome(int start, int end, byte[] template) {
    mStart = start;
    mEnd = end;
    mTemplate = template;
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
    final int ix = index + mStart;
    if (ix < 0 || ix >= mTemplate.length) {
      return 0;
    }
    return mTemplate[ix];
  }

  @Override
  public int subsequenceLength() {
    return mEnd - mStart;
  }

  @Override
  public int templateLength() {
    return mTemplate.length;
  }

  @Override
  public boolean isInverted() {
    return false;
  }

}
