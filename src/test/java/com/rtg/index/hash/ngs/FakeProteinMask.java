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
package com.rtg.index.hash.ngs;

import java.io.IOException;

import org.junit.Assert;

import com.rtg.index.hash.ngs.general.Skeleton;
import com.rtg.index.hash.ngs.protein.ProteinMask;
import com.rtg.launcher.HashingRegion;

class FakeProteinMask extends ProteinMask {

  int mReadCalls;
  int mTemplateCalls;
  int mClones;

  FakeProteinMask(Skeleton sk, ReadCall readCall, TemplateCall templateCall) {
    super(sk, readCall, templateCall);
  }

  @Override
  public void readAll(int readId, boolean reverse) {
    Assert.assertFalse(reverse);
    ++mReadCalls;
  }

  @Override
  public void templateForward(int endPosition) {
    ++mTemplateCalls;
  }

  @Override
  public NgsHashFunction threadClone(HashingRegion region) throws IOException {
    ++mClones;
    return super.threadClone(region);
  }
}
