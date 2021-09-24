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
