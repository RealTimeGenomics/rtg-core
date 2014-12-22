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

/**
 * Writes all results to specified <code>Appendable</code>.
 *
 */
public class RawPositionWriter implements PositionWriter {

  private static class SharedRawVars {
    Appendable mOut;
    Appendable mUnmapped;
    double mTotalScore = 0.0;
    boolean mSeen = false;
    boolean mHeaderWritten = false;
    SurrogateRegion mSurrogate = null;
  }

  private final SharedRawVars mSharedVars;
  private final boolean mSuppressEndQuery;

  /**
   * Constructor
   * @param out stream to write results to
   * @param unmapped stream to write unmapped to
   */
  public RawPositionWriter(final Appendable out, final Appendable unmapped) {
    mSharedVars = new SharedRawVars();
    mSharedVars.mOut = out;
    mSharedVars.mUnmapped = unmapped;
    mSuppressEndQuery = false;
  }


  @Override
  public void endQuery(final int queryId) throws IOException {
    if (!mSuppressEndQuery) {
      if (!mSharedVars.mSeen && mSharedVars.mUnmapped != null) {
        mSharedVars.mUnmapped.append(String.valueOf(queryId));
        mSharedVars.mUnmapped.append(com.rtg.util.StringUtils.LS);
      }
      mSharedVars.mSeen = false;
    }
  }

  @Override
  public void write(final AbstractGappedRegion<?> region, final int sequenceId, final Frame queryFrame, final int queryLength, final int queryEffectiveLength) throws IOException {
    //System.err.println("queryLength=" + queryLength + " queryEffectiveLength=" + queryEffectiveLength + " RawPositionWriter");
    mSharedVars.mSurrogate = region.surrogateOutput(mSharedVars.mSurrogate, sequenceId, queryFrame, queryLength, queryEffectiveLength);
    if (!mSharedVars.mHeaderWritten) {
      mSharedVars.mHeaderWritten = true;
      mSharedVars.mSurrogate.writeHeader(mSharedVars.mOut);
    }
    if (mSharedVars.mOut != null) {
      mSharedVars.mSeen |= mSharedVars.mSurrogate.write(mSharedVars.mOut, null, null);
    }
    mSharedVars.mTotalScore += mSharedVars.mSurrogate.score();
  }

  @Override
  public double score() {
    return mSharedVars.mTotalScore;
  }

  @Override
  public void endAll() {
    //do nothing
  }

  /**
   * Return a clone of this suitable for use in the reverse frame
   * @return the clone
   */
  @Override
  public PositionWriter reverseClone() {
    throw new UnsupportedOperationException();
  }


}
