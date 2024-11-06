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

package com.rtg.assembler.graph.implementation;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.PathsIterator;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
@TestClass({"com.rtg.assembler.graph.implementation.GraphImplementationTest", "com.rtg.assembler.graph.implementation.IteratorLocalTest"})
public final class IteratorLocal extends IntegralAbstract implements PathsIterator {

  private boolean mShowDeleted;

  private long mContigId;

  private long mNext = -1;

  private int mIndex = -1;

  private final GraphImplementation mGraphImplementation;

  /**
   * @param graphImplementation the graph that is parent of this (was an inner class before being moved out)
   */
  public IteratorLocal(GraphImplementation graphImplementation) {
    mGraphImplementation = graphImplementation;
  }

  @Override
  public boolean integrity() {
    if (mNext == -1) {
      return true;
    }
    mGraphImplementation.absContig(mContigId);
    Exam.assertTrue(-1 <= mIndex);
    return true;
  }

  @Override
  public Graph graph() {
    return mGraphImplementation;
  }

  @Override
  public void set(long contigId) {
    set(contigId, false);
  }

  @Override
  public void set(long contigId, boolean showDeleted) {
    mContigId = contigId;
    final long acontig = mGraphImplementation.absContig(contigId);
    mNext = mGraphImplementation.mLastPath.get(acontig);
    mIndex = -1;
    mShowDeleted = showDeleted;
  }

  @Override
  public int contigIndex() {
    if (mIndex == -1) {
      throw new IllegalStateException();
    }
    return mIndex;
  }

  @Override
  public long nextPathId() {
    while (true) {
      if (mNext == -1) {
        mIndex = -1;
        return 0;
      }
      final long curr = mNext;
      mNext = mGraphImplementation.mPathIndexPrev.get(curr);
      final long offset = mGraphImplementation.mPathIndexOffset.get(curr);
      final long path = mGraphImplementation.mPathIndex.get(curr - offset);
      if (!mShowDeleted && mGraphImplementation.pathDeleted(path)) {
        continue;
      }
      final long pcontig = mGraphImplementation.mPathIndex.get(curr);
      final long spath;
      if (pcontig != mContigId) {
        assert pcontig == -mContigId;
        mIndex  = mGraphImplementation.pathLength(path) - (int) offset;
        spath = -path;
      } else {
        mIndex  = (int) offset - 1;
        spath = path;
      }
      return spath;
    }
  }
}
