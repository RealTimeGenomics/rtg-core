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
