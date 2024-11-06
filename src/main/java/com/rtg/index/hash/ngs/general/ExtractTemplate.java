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
package com.rtg.index.hash.ngs.general;

import java.io.IOException;

import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.util.integrity.Exam;

/**
 */
public final class ExtractTemplate extends ExtractAbstract {

  final TemplateCall mTemplateCall;

  private final int mIndex;

  private int mEndPosition = -1;

  ExtractTemplate(final SingleMask singleMask, final TemplateCall templateCall, final int index) {
    super(singleMask);
    mTemplateCall = templateCall;
    mIndex = index;
    integrity();
  }

  void templateCall(final int endPosition, final long v0, final long v1) throws IOException {
    mEndPosition = endPosition;
    //this ends up calling masked() below and uses the mReadId
    maskIndel(v0, v1);
  }

  @Override
  protected void masked(final long bits) throws IOException {
    //System.err.println("hash=" + com.rtg.util.Utils.toBitsSep(bits) + " " + bits);
    mTemplateCall.templateCall(mEndPosition, bits, mIndex);
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("ExtractTemplate templateCall=").append(mTemplateCall).append(" index=").append(mIndex).append(" end=").append(mEndPosition);
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mTemplateCall != null);
    Exam.assertTrue(mIndex >= 0);
    return true;
  }
}
