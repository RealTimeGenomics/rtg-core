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
package com.rtg.index.hash.ngs.protein;

import java.io.IOException;

import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.general.SingleMask;
import com.rtg.util.integrity.Exam;

/**
 */
public final class ProteinExtractRead extends AbstractProteinExtract {

  private final ReadCall mReadCall;

  private final int mIndex;

  private int mReadId = -1;

  ProteinExtractRead(final SingleMask singleMask, final ReadCall readCall, final int index) {
    super(singleMask);
    mReadCall = readCall;
    mIndex = index;
    integrity();
  }

  void readCall(int id, long[] vs) throws IOException {
    mReadId = id;
    //this ends up calling masked() below and uses the mReadId
    mask(vs);
  }

  @Override
  protected void masked(final long bits) {
    if (bits == 0) { //if resulting hash is zero do not add into index, in protein in translate to all X's
      return;
    }
    mReadCall.readCall(mReadId, bits, mIndex);
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("ProteinExtractRead readCall=").append(mReadCall).append(" index=").append(mIndex).append(" id=").append(mReadId);
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mReadCall != null);
    Exam.assertTrue(mIndex >= 0);
    return true;
  }
}
