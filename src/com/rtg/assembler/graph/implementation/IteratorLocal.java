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
