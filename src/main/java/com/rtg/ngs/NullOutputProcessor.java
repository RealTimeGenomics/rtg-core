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
package com.rtg.ngs;

import java.util.ArrayList;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Does absolutely nothing but counting of hits, to allow profiling of index searching.
 *
 */
public class NullOutputProcessor implements OutputProcessor {

  private long mHitCount = 0;
  private final ArrayList<NullOutputProcessor> mChildren = new ArrayList<>();

  /**
   * Does nothing
   */
  public NullOutputProcessor() {
  }

  @Override
  public void process(final long templateId, final String frame, final int readId, final int tStart, final int score, final int scoreIndel) {
    ++mHitCount;
    /*
    System.out.println(""
                + templateId + "\t"
                + frame + "\t"
                + readId + "\t"
                + tStart + "\t"
                + score + "\t"
                + scoreIndel
                );
     */
  }

  @Override
  public void finish() {
    long totalHits = mHitCount; // Single thread case?
    for (NullOutputProcessor o : mChildren) {
      totalHits += o.mHitCount;
    }
    Diagnostic.developerLog("NullOutputProcessor total number of hits:" + totalHits);
  }
  @Override
  public OutputProcessor threadClone(HashingRegion region) {
    final NullOutputProcessor inner = new NullOutputProcessor();
    mChildren.add(inner);
    Diagnostic.developerLog("NullOutputProcessor created for region:" + region);
    return inner;
  }

  @Override
  public void threadFinish() {
  }
  @Override
  public void close() {
  }

  @Override
  public String toString() {
    return "NullOutputProcessor";
  }
}

