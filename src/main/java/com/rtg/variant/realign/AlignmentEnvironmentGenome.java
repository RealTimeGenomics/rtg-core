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
