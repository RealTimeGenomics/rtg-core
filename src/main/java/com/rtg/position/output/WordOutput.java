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
