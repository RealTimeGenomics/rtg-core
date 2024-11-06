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

import com.rtg.index.Finder;
import com.rtg.index.IndexSet;
import com.rtg.launcher.HashingRegion;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Computes minimum distance between candidate probes and template sequences..
 */
public class PcrTemplateCallImplementation extends IntegralAbstract implements TemplateCall {

  private NgsHashFunction mHashFunction;

  private final IndexSet mIndexes;

  private final byte[] mDistances;

  private Finder mHit;

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("PcrTemplateCallImplementation");
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mIndexes != null && mIndexes.size() > 0);
    Exam.assertTrue(mDistances != null && mDistances.length > 0);
    Exam.assertTrue(mHit != null);
    return true;
  }

  /**
   * @param indexes indexes selected by the hash function.
   * @param distances minimum match distance so far
   */
  public PcrTemplateCallImplementation(final IndexSet indexes, final byte[] distances) {
    mHit = makeHit();
    mDistances = distances;
    mIndexes = indexes;
    //System.err.println("numberReads=" + mNumberReads + " threshold=" + threshold);
  }


  private Finder makeHit() {
    return new Finder() {
      @Override
      public boolean found(final long readIdl) {
        final int readId = (int) readIdl;
        final int score = mHashFunction.fastScore(readId);
        //System.err.println(" readId=" + readId + " score=" + score);
        assert score > 0 && score <= Byte.MAX_VALUE : score;
        if (score < mDistances[readId]) {
          mDistances[readId] = (byte) score;
        }
        return true;
      }
    };
  }

  @Override
  public PcrTemplateCallImplementation clone() throws CloneNotSupportedException {
    final PcrTemplateCallImplementation clone = (PcrTemplateCallImplementation) super.clone();
    clone.mHit = clone.makeHit();
    assert integrity();
    return clone;
  }

  @Override
  public TemplateCall threadClone(final HashingRegion region) {
    if (region != HashingRegion.NONE) {
      throw new UnsupportedOperationException();
    }
    try {
    return clone();
    } catch (final CloneNotSupportedException e) {
      throw new RuntimeException(e);
    }
  }

  @Override
  public void threadFinish() {
  }

  @Override
  public void set(final long name, final int length) {
    //do nothing
  }

  @Override
  public void setReverse(final boolean reverse) {
    // do nothing
  }

  @Override
  public boolean isReverse() {
    return false;
  }

  @Override
  public void setHashFunction(final NgsHashFunction hashFunction) {
    mHashFunction = hashFunction;
  }

  @Override
  public void done() {
    //do nothing
  }

  @Override
  public void endSequence() {
    //do nothing
  }

  @Override
  public void templateCall(final int endPosition, final long hash, final int index) throws IOException {
    //System.err.println("search index=" + index + " hash=" + Utils.toBits(hash));
    mIndexes.get(index).search(hash, mHit);
  }

  @Override
  public void logStatistics() {
    // do nothing
  }
}
