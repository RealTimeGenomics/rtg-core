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
  private String mTemplateName;


  CalibratorCigarParser(Calibrator calibrator) {
    mCalibrator = calibrator;
  }
  void setTemplateName(String name) {
    mTemplateName = name;
  }
  boolean include() {
    return mCalibrator.inRange(mTemplateName, getTemplatePosition());
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
   *  Not thread safe (due to <code>getOrCreateHistogram</code> calls)
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
      mMismatchCount++;
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

  @Override
  protected void doUnknownOnTemplate(int readNt, int templateNt) throws BadSuperCigarException {
    // Skip sites where the template is N, we can then use Ns to mask sites where calibration should
    // not be carried out (e.g. at sites of known variants)
  }

  @Override
  protected void doUnknownOnRead() throws BadSuperCigarException {
    // N's on the read are simply skipped
  }
}
