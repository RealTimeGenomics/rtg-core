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

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;

import com.rtg.mode.Frame;
import com.rtg.reader.NamesInterface;

/**
 * Surrogate for <code>GappedRegion</code>.
 *
 */
class SurrogateGappedRegion implements SurrogateRegion {

  private final Frame[] mSubjectFrames;
  private final int mNumberFrames;

  private int mQueryId;
  private Frame mQueryFrame;
  private int mQueryStart;
  private int mQueryEnd;
  private int mBuildSeqId;
  private int mBuildStart;
  private int mBuildEnd;

  final SurrogateRegion initialize(final AbstractGappedRegion<?> region, final int queryId, final Frame queryFrame) {
    mQueryId = queryId;
    mQueryFrame = queryFrame;
    mQueryStart = region.mQueryStart;
    mQueryEnd = region.queryEnd();
    mBuildSeqId = region.sequenceId();
    mBuildStart = region.mBuildStart;
    mBuildEnd = region.mBuildEnd;
    return this;
  }

  SurrogateGappedRegion(final AbstractGappedRegion<?> region, final int queryId, final Frame queryFrame, final Frame[] subjectFrames) {
    mSubjectFrames = subjectFrames;
    mNumberFrames = mSubjectFrames.length;
    initialize(region, queryId, queryFrame);
  }

  @Override
  public boolean write(final Appendable out, final NamesInterface subjectNames, final NamesInterface queryNames) throws IOException {
    final Frame subjectFrame = mSubjectFrames[mBuildSeqId % mNumberFrames];
    final long sequenceId = mBuildSeqId / mNumberFrames;
    final String queryId = queryNames == null ? Long.toString(queryId()) : queryNames.name(queryId());
    final String subjectId = subjectNames == null ? Long.toString(sequenceId) : subjectNames.name(sequenceId);
      out.append("").append(queryId).append("\t").append(mQueryFrame.display()).append("\t").append(String.valueOf(mQueryStart + 1)).append("\t").append(String.valueOf(mQueryEnd + 1)).append("\t").append(subjectId).append("\t").append(subjectFrame.display()).append("\t").append(String.valueOf(mBuildStart + 1)).append("\t").append(String.valueOf(mBuildEnd + 1)).append(LS);
    return true;
  }

  @Override
  public boolean scoreAllowed() {
    return true;
  }



  /**
   * @see com.rtg.position.output.SurrogateRegion#writeHeader(java.lang.Appendable)
   */
  @Override
  public void writeHeader(final Appendable out) throws IOException {
    out.append("#" + "query-id\t" + "query-frame\t" + "query-start\t" + "query-end\t" + "subject-id\t" + "subject-frame\t" + "subject-start\t" + "subject-end").append(LS);
  }

  @Override
  public double score() {
    return mBuildEnd - mBuildStart + 1;
  }

  @Override
  public int queryId() {
    return mQueryId;
  }
}
