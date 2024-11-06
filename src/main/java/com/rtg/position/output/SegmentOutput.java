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
package com.rtg.position.output;

import java.io.IOException;

import com.rtg.launcher.BuildParams;
import com.rtg.mode.Frame;
import com.rtg.util.integrity.Exam;

/**
 */
public class SegmentOutput extends AbstractPositionOutput implements SegmentWriter {

  private int mIndexPosition = NULL;

  protected final int mWindowSize;

  protected final int mStepSize;

  private final int mThreshold;

  private final SegmentCollector mCollector;

  private final SegmentCollection mFree;

  private final SegmentCollection[] mCollections;

  private final int[] mPositions;

  private SharedSegmentVars mSharedVars;

  private SegmentCollection mCurrentCollection;

  private static class SharedSegmentVars {
    boolean mHeaderWritten = false;

    long mScore;
  }

  /**
   * @param params parameters.
   * @param out where to write the output.
   */
  public SegmentOutput(final PositionParams params, final Appendable out) {
    super(params, out);
    final BuildParams buildParams = params.build();
    mSharedVars = new SharedSegmentVars();
    mWindowSize = buildParams.windowSize();
    mStepSize = buildParams.stepSize();
    mThreshold = params.hashCountThreshold();
    mCollector = new SegmentCollector(mThreshold, buildParams.windowSize(), mStepSize, this);
    mCurrentCollection = new SegmentCollection(params.hashCountThreshold());
    mCollections = new SegmentCollection[mStepSize];
    mPositions = new int[mStepSize];
    for (int i = 0; i < mCollections.length; ++i) {
      mCollections[i] = new SegmentCollection(params.hashCountThreshold());
      mPositions[i] = NULL;
    }
    final int free = mThreshold * (mStepSize + 1);
    mFree = new SegmentCollection(free);
    for (int i = 0; i < free; ++i) {
      mFree.add(new Segment());
    }
    integrity();
    assert globalIntegrity();
  }

  private SegmentOutput(final PositionParams params, final Appendable out, final SharedSegmentVars sharedVars) {
    this(params, out);
    mSharedVars = sharedVars;
  }



  @Override
  public void hit(final int seqId, final int posn) throws IOException {
    final int p = posn + mWindowSize - 1;
    super.hit(seqId, p);
    mCollector.add(seqId, p);
  }

  @Override
  public void setPosition(final int position) throws IOException {
    //System.err.println("SO.setPosition position=" + position);
    if (mSearchPosition >= mStepSize) {
      final int delta =  position - mSearchPosition;
      //assert delta > 0 : mSearchPosition + ":" + position;
      final int end = delta > mStepSize ? mSearchPosition + mStepSize + 1 : position;
      //System.err.println("SO.setPosition mSearchPosition=" + mSearchPosition + " delta=" + delta + " end=" + end);
      for (int i = mSearchPosition + 1, index = i % mStepSize; i < end; ++i, ++index) {
        if (index == mStepSize) {
          index = 0;
        }
        mCollections[index].flush(this, mFree, mPositions[index]);
        mPositions[index] = NULL;
      }
    }
    mIndexPosition = position % mStepSize;
    super.setPosition(position);
  }

  @Override
  public void endPosition() throws IOException {
    //System.err.println("SO.endPosition");
    assert mIndexPosition >= 0 && mIndexPosition < mStepSize;
    assert mSearchPosition >= 0;
    super.endPosition();
    //process all the accumulated hits
    mCollector.endPosition(mCollections[mIndexPosition], mCurrentCollection, mFree, mPositions[mIndexPosition]);
    final SegmentCollection temp = mCollections[mIndexPosition];
    mCollections[mIndexPosition] = mCurrentCollection;
    mCurrentCollection = temp;
    mPositions[mIndexPosition] = mCollections[mIndexPosition].size() == 0 ? NULL : mSearchPosition;
    assert mCurrentCollection.size() == 0;
    assert globalIntegrity();
  }

  @Override
  public void endQuery() throws IOException {
    //System.err.println("SO.endQuery");
    //try and keep the outputs in order
    for (int i = mSearchPosition + 1, index = i % mStepSize; i < mSearchPosition + mStepSize + 1; ++i, ++index) {
      if (index == mStepSize) {
        index = 0;
      }
      //System.err.println("SO.endQuery.flush i=" + i + " index=" + index);
      mCollections[index].flush(this, mFree, mPositions[index]);
      mPositions[index] = NULL;
    }

    super.endQuery();
    mIndexPosition = NULL;
    assert globalIntegrity();
  }

  @Override
  public void write(final Segment segment, final int searchPosition) throws IOException {
    final int seqId = segment.seqId();
    final int posn = segment.start();
    final int length = segment.end() - posn + 1;
    mSharedVars.mScore += length;
    final long sequenceId = seqId / mNumberSFrames;
    final Frame sFrame = mSubjectFrames[seqId % mNumberSFrames];
    if (!mSharedVars.mHeaderWritten) {
      mSharedVars.mHeaderWritten = true;
      writeHeader();
    }
    mOut.append("").append(String.valueOf(mSearchSeqId)).append("\t").append(mSearchFrame.display()).append("\t").append(String.valueOf(searchPosition - length + mWindowSize + 1)).append("\t").append(String.valueOf(sequenceId)).append("\t").append(sFrame.display()).append("\t").append(String.valueOf(posn + 1)).append("\t").append(String.valueOf(length)).append(com.rtg.util.StringUtils.LS);
  }

  private void writeHeader() throws IOException {
    mOut.append("#query-id\t" + "query-frame\t" + "query-start\t" + "subject-id\t" + "subject-frame\t" + "subject-start\t" + "length").append(com.rtg.util.StringUtils.LS);

  }

  @Override
  public double score() {
    return mSharedVars.mScore;
  }

  @Override
  public long bytes() {
    return 0;
  }

  @Override
  public void memToString(final StringBuilder sb) {
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("SegmentOutput");
  }

  @Override
  public final boolean globalIntegrity() {
    integrity();
    int sum = 0;
    for (int i = 0; i < mCollections.length; ++i) {
      final SegmentCollection coll = mCollections[i];
      coll.globalIntegrity();
      sum += coll.size();
      final int posn = mPositions[i];
      Exam.assertTrue(posn >= NULL);
      if (posn == NULL) {
        Exam.assertEquals(0, coll.size());
      }
    }
    mFree.globalIntegrity();
    sum += mFree.size();
    Exam.assertEquals(mThreshold * (mStepSize + 1), sum);
    mCollector.globalIntegrity();
    if (mSearchSeqId == NULL) {
      for (final SegmentCollection coll : mCollections) {
        coll.globalIntegrity();
        Exam.assertEquals(0, coll.size());
      }
    }
    return true;
  }

  @Override
  public final boolean integrity() {
    super.integrity();
    Exam.assertTrue(mSubjectFrames != null && mNumberSFrames == mSubjectFrames.length && mNumberSFrames > 0);
    Exam.assertTrue(mOut != null);
    Exam.assertTrue(mStepSize > 0);
    Exam.assertTrue(mThreshold > 0);
    if (mSearchSeqId == NULL || mSearchFrame == null) {
      Exam.assertEquals(NULL, mIndexPosition);
      for (int i = 0; i < mStepSize; ++i) {
        Exam.assertEquals(NULL, mPositions[i]);
        Exam.assertEquals(0, mCollections[i].size());
      }
    } else {
      Exam.assertTrue(mIndexPosition >= 0 && mIndexPosition < mStepSize);
      Exam.assertTrue(mSearchPosition == NULL || mSearchPosition % mStepSize == mIndexPosition);
      for (int i = 0; i < mStepSize; ++i) {
        Exam.assertEquals(NULL == mPositions[i], 0 == mCollections[i].size());
      }
    }
    Exam.assertTrue(mCollector != null);
    Exam.assertTrue(mFree != null);
    if (mPositions == null) {
      throw new NullPointerException();
    } else {
      Exam.assertTrue(mStepSize == mPositions.length);
      if (mCollections == null) {
        throw new NullPointerException();
      } else {
        Exam.assertTrue(mStepSize == mCollections.length);
        for (int i = 0; i < mCollections.length; ++i) {
          Exam.assertTrue(mCollections[i] != null);
          Exam.assertTrue(mPositions[i] >= NULL);
        }
      }
      Exam.assertTrue(mCurrentCollection != null && 0 == mCurrentCollection.size());
      Exam.assertTrue(mSharedVars.mScore >= 0);
      return true;
    }
  }

  @Override
  public PositionOutput reverseClone() {
    return new SegmentOutput(mParams, mOut, mSharedVars);
  }


}
