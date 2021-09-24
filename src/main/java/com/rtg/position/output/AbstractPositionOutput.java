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
package com.rtg.position.output;

import java.io.IOException;

import com.rtg.launcher.BuildParams;
import com.rtg.mode.Frame;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * The calls obey the following grammar:
 * <br>
 * <code> (nextSeq (nextQuery (setPosition hit* endPosition)* endQuery)* endQuerySequence) endAll</code>
 */
public abstract class AbstractPositionOutput extends IntegralAbstract implements PositionOutput {

  protected static final int NULL = -1;

  protected enum State {
    SEQUENCE, QUERY, POSITION, HIT
  }

  protected final PositionParams mParams;

  protected final Appendable mOut;

  protected final Frame[] mSubjectFrames;

  protected final int mNumberSFrames;

  private final int mWordSize;

  /** Number of residues in the query sequence. */
  protected int mQueryLength = NULL;

  /** Number of windows in the query sequence. Used for scoring. */
  protected int mQueryEffectiveLength = NULL;

  protected Frame mSearchFrame = null;

  protected int mSearchSeqId = NULL;

  protected int mSearchPosition = NULL;

  protected boolean mQuerySequence = false;

  protected State mState = State.SEQUENCE;

  /**
   * @param params parameters.
   * @param out where to write the output.
   */
  protected AbstractPositionOutput(final PositionParams params, final Appendable out) {
    mParams = params;
    final BuildParams buildParams = params.build();
    final Frame[] subjectFrames = buildParams.sequences().mode().allFrames();
    mSubjectFrames = subjectFrames.clone();
    mNumberSFrames = subjectFrames.length;
    mWordSize = buildParams.windowSize();
    mOut = out;
  }

  /**
   * This version for testing only.
   * @param subjectFrames frames used for subjects.
   * @param out where to write the output.
   */
  protected AbstractPositionOutput(final Frame[] subjectFrames, final Appendable out) {
    mParams = null;
    mSubjectFrames = subjectFrames.clone();
    mNumberSFrames = subjectFrames.length;
    mWordSize = 1; //doesnt matter for the few tests that use this
    mOut = out;
  }

  @Override
  public void hit(final int seqId, final int posn) throws IOException {
    assert mState == State.HIT;
  }

  @Override
  public void setPosition(final int position)  throws IOException {
    assert mState == State.POSITION;
    mSearchPosition = position;
    mState = State.HIT;
  }

  @Override
  public void endPosition() throws IOException {
    assert mState == State.HIT;
    mState = State.POSITION;
  }

  @Override
  public void nextSequence(int seqId, final int length) throws IOException {
    assert mState == State.SEQUENCE;
    mState = State.QUERY;
    mQueryLength = length;
    final int l = length - mWordSize + 1;
    mQueryEffectiveLength = l < 0 ? 0 : l;
    mSearchSeqId = seqId;
  }

  @Override
  public void nextSequence(int seqId, final int length, final int usedLength, final byte[] sequence) throws IOException {
    assert length >= usedLength && usedLength >= 0;
    nextSequence(seqId, length);
  }

  @Override
  public void nextQuery(final Frame frame, final int seqId) {
    assert mState == State.QUERY;
    mState = State.POSITION;
    mSearchFrame = frame;
    mSearchSeqId = seqId;
    mQuerySequence = true;
  }

  @Override
  public void endQuery() throws IOException {
    assert mState == State.POSITION;
    mSearchPosition = NULL;
    mSearchFrame = null;
    mState = State.QUERY;
  }

  @Override
  public void endQuerySequence() throws IOException {
    assert mState == State.QUERY;
    mQuerySequence = false;
    mSearchSeqId = NULL;
    mQueryLength = NULL;
    mQueryEffectiveLength = NULL;
    mState = State.SEQUENCE;
  }

  @Override
  public void endAll() throws IOException {
    assert mState == State.SEQUENCE;
  }

  /**
   * Create a human readable string that analyses the memory usage
   * of the current data structure.
   * @param sb the string describing memory usage is placed here.
   */
  public abstract void memToString(StringBuilder sb);

  /**
   * Create a human readable string that analyses the memory usage
   * of the current data structure.
   * @return the string describing memory usage.
   */
  public String memToString() {
    final StringBuilder sb = new StringBuilder();
    memToString(sb);
    return sb.toString();
  }

  /**
   * Compute the total number of bytes in the current data structure.
   * May ignore terms of O(1).
   * @return the total number of bytes in the current data structure.
   */
  public abstract long bytes();

  /**
   * Create a human readable string that analyses the memory usage
   * of the current data structure.
   * @param sb the string describing memory usage is placed here.
   */
  @Override
  public abstract void toString(StringBuilder sb);

  @Override
  public boolean integrity() {
    Exam.assertTrue(mWordSize >= 1);
    Exam.assertTrue(mState != null);
    Exam.assertTrue(mSubjectFrames != null && mNumberSFrames == mSubjectFrames.length && mNumberSFrames > 0);

    if (mState == State.SEQUENCE) {
      Exam.assertEquals(NULL, mQueryLength);
      Exam.assertEquals(NULL, mQueryEffectiveLength);
      Exam.assertEquals(NULL, mSearchPosition);
      Exam.assertTrue(!mQuerySequence);
      Exam.assertEquals(NULL, mSearchSeqId);
      Exam.assertEquals(null, mSearchFrame);
    } else if (mState == State.QUERY) {
      Exam.assertEquals(NULL, mSearchPosition);
      Exam.assertTrue(mQueryLength >= 0); //allow 0 for some tests
      Exam.assertTrue(mQueryEffectiveLength >= 0);
      if (mQuerySequence) {
        Exam.assertTrue(mSearchSeqId >= 0);
      } else {
        Exam.assertEquals(NULL, mSearchSeqId);
      }
      Exam.assertEquals(null, mSearchFrame);
    } else {
      Exam.assertTrue(mSearchPosition >= 0);
      Exam.assertTrue(mSearchSeqId >= 0);
      Exam.assertTrue(mSearchFrame != null);
    }
    return true;
  }

}
