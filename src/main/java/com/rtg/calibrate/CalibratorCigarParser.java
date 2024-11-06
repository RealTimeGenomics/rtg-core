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

package com.rtg.calibrate;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.mode.DnaUtils;
import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.SuperCigarParser;
import com.rtg.util.Histogram;

/**
 * Updates calibration counts during CIGAR parsing.
 *
 */
@TestClass("com.rtg.calibrate.CalibratorTest")
class CalibratorCigarParser extends SuperCigarParser {

  private final Calibrator mCalibrator;
  private String mCurrentReadGroup = "";
  private Histogram mInsHist;
  private Histogram mDelHist;
  private Histogram mCGOverHist;
  private Histogram mCGGapHist;
  private Histogram mMnpHist;


  CalibratorCigarParser(Calibrator calibrator) {
    mCalibrator = calibrator;
  }

  boolean include() {
    return mCalibrator.inRange(getTemplatePosition());
  }

  private static final int DEFAULT_QUALITY = 20;

  private int mMismatchCount = 0;

  private boolean mGotohAligned = false;
  private byte[] mQualities = null;

  private void setKeys() {
    // NOTE value of mReadGroup in mCalibrator is NOT immutable, hence need to reset these
    final String readGroup = mCalibrator.mReadGroup;
    if (!mCurrentReadGroup.equals(readGroup)) {
      mCurrentReadGroup = readGroup == null ? "" : readGroup;
      mInsHist = null;
      mDelHist = null;
      mCGOverHist = null;
      mCGGapHist = null;
      mMnpHist = null;
    }
  }

  @Override
  public void setCigar(String superCigar, String readDelta) {
    mGotohAligned = true;
    super.setCigar(superCigar, readDelta);
    mQualities = null;
    setKeys();
  }

  @Override
  public void setStandardCigar(String cigar, byte[] read, int readLength) {
    mGotohAligned = false;
    super.setStandardCigar(cigar, read, readLength);
    mQualities = null;
    setKeys();
  }

  /**
   * Set the expanded qualities for this SAM record
   * @param qualities the expanded qualities
   */
  void setQualities(byte[] qualities) {
    mQualities = qualities;
  }

  int getCurrentQuality() throws BadSuperCigarException {
    if (mQualities == null) {
      return DEFAULT_QUALITY;
    }
    if (getReadPosition() >= mQualities.length) {
      throw new BadSuperCigarException("readPos " + getReadPosition() + " > qual.len in SAM record");
    }
    return mQualities[getReadPosition()];
  }

  /**
   * Not thread safe (due to <code>getOrCreateHistogram</code> calls)
   * @throws BadSuperCigarException on bad cigar
   */
  @Override
  protected void doChunk(char ch, int count) throws BadSuperCigarException {
    startMatchMismatch();
    super.doChunk(ch, count);
    if (include()) {
      switch (ch) {
        case 'I':
          if (mInsHist == null) {
            mInsHist = mCalibrator.getOrCreateHistogram(Calibrator.INS_DIST + ":" + mCurrentReadGroup);
          }
          mInsHist.increment(count);
          break;
        case 'D':
          if (mDelHist == null) {
            mDelHist = mCalibrator.getOrCreateHistogram(Calibrator.DEL_DIST + ":" + mCurrentReadGroup);
          }
          mDelHist.increment(count);
          break;
        case 'B':
          if (mGotohAligned) {
            if (mCGOverHist == null) {
              mCGOverHist = mCalibrator.getOrCreateHistogram(Calibrator.CGOVER_DIST + ":" + mCurrentReadGroup);
            }
            mCGOverHist.increment(count);
          }
          break;
        case 'N':
          if (mGotohAligned) {
            if (mCGGapHist == null) {
              mCGGapHist = mCalibrator.getOrCreateHistogram(Calibrator.CGGAP_DIST + ":" + mCurrentReadGroup);
            }
            mCGGapHist.increment(count);
          }
          break;
        default:
          break;
      }
    }
    endMatchMismatch();
  }

  private void endMatchMismatch() {
    if (include()) {
      if (mMismatchCount != 0) {
        if (mMnpHist == null) {
          mMnpHist = mCalibrator.getOrCreateHistogram(Calibrator.MNP_DIST + ":" + mCurrentReadGroup);
        }
        mMnpHist.increment(mMismatchCount);
      }
      mMismatchCount = 0;
    }
  }

  private void startMatchMismatch() {
    if (include()) {
      mMismatchCount = 0;
    }
  }

  @Override
  protected void doReadOnly(int readNt) throws BadSuperCigarException {
    if (include()) {
      if (readNt != DnaUtils.UNKNOWN_RESIDUE) {
        mCalibrator.findStats(this).seenInsert(1);
      }
    }
  }

  @Override
  protected void doTemplateOnly(int templateNt) throws BadSuperCigarException {
    if (include()) {
      mCalibrator.findStats(this).seenDelete(1);
    }
  }

  @Override
  protected void doSubstitution(int readNt, int templateNt) throws BadSuperCigarException {
    if (include()) {
      if (readNt != DnaUtils.UNKNOWN_RESIDUE) {
        mCalibrator.findStats(this).seenMnp(1);
      }
      ++mMismatchCount;
    }
  }

  @Override
  protected void doEquality(int readNt, int nt) throws BadSuperCigarException {
    if (include()) {
      if (readNt != DnaUtils.UNKNOWN_RESIDUE) {
        mCalibrator.findStats(this).seenEquals(1);
      }
      endMatchMismatch();
      startMatchMismatch();
    }
  }
}
