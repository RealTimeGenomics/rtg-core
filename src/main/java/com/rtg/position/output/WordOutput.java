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

import com.rtg.mode.Frame;
import com.rtg.util.integrity.Exam;

/**
 */
public class WordOutput extends AbstractPositionOutput {

  private static class SharedWordOutputVars {
    Appendable mUnmapped;
    boolean mHeaderWritten = false;
    long mScore;
    boolean mSeen = false;
  }

  private final boolean mSuppressEndQuery;

  private final SharedWordOutputVars mSharedVars;

  /**
   * This version for testing only.
   * @param subjectFrames frames used for subjects.
   * @param out where to write the output.
   * @param unmappedOut destination for unmapped sequences
   */
  WordOutput(final Frame[] subjectFrames, final Appendable out, final Appendable unmappedOut) {
    super(subjectFrames, out);
    mSharedVars = new SharedWordOutputVars();
    mSharedVars.mUnmapped = unmappedOut;
    integrity();
    mSuppressEndQuery = false;
  }

  /**
   * Constructor for <code>reverseClone</code>
   * @param params as in above
   * @param out as above
   * @param sharedVars the shared variables
   */
  private WordOutput(final PositionParams params, final Appendable out, final SharedWordOutputVars sharedVars, final boolean suppressEndQuery) {
    super(params, out);
    this.mSharedVars = sharedVars;
    mSuppressEndQuery = suppressEndQuery;
  }

  @Override
  public void hit(final int seqId, final int posn) throws IOException {
    super.hit(seqId, posn);
    mSharedVars.mScore++;
    final long sequenceId = seqId / mNumberSFrames;
    final Frame sFrame = mSubjectFrames[seqId % mNumberSFrames];
    mSharedVars.mSeen = true;
    if (!mSharedVars.mHeaderWritten) {
      mSharedVars.mHeaderWritten = true;
      writeHeader();
    }
    mOut.append("").append(String.valueOf(mSearchSeqId)).append("\t").append(mSearchFrame.display()).append("\t").append(String.valueOf(mSearchPosition + 1)).append("\t").append(String.valueOf(sequenceId)).append("\t").append(sFrame.display()).append("\t").append(String.valueOf(posn + 1)).append(com.rtg.util.StringUtils.LS);
  }

  @Override
  public void endQuerySequence() throws IOException {
    if (!mSuppressEndQuery) {
      if (!mSharedVars.mSeen && mSharedVars.mUnmapped != null) {
        mSharedVars.mUnmapped.append(String.valueOf(mSearchSeqId)).append(com.rtg.util.StringUtils.LS);
      }
      mSharedVars.mSeen = false;
    }
    super.endQuerySequence();
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
    sb.append("WordOutput");
  }

  private void writeHeader() throws IOException {
    mOut.append("#" + "query-id\t" + "query-frame\t" + "query-start\t" + "subject-id\t" + "subject-frame\t" + "subject-start").append(com.rtg.util.StringUtils.LS);
  }

  @Override
  public final boolean integrity() {
    super.integrity();
    Exam.assertTrue(mSharedVars.mUnmapped != null);
    Exam.assertTrue(mSharedVars.mScore >= 0);
    return true;
  }

  @Override
  public PositionOutput reverseClone() {
    return new WordOutput(mParams, mOut, mSharedVars, true);
  }

}
