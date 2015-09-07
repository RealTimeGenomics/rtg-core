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
package com.rtg.sam;

import java.util.Arrays;

import com.rtg.alignment.CgGotohEditDistance;
import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.CgUtils;
import com.rtg.reader.FastaUtils;
import com.rtg.util.Utils;

import htsjdk.samtools.SAMRecord;


/**
 * Validates a super cigar against the template and original read
 */
public class SuperCigarValidator extends SuperCigarParser {

  private SAMRecord mSamRecord = null;
  private byte[] mSdfRead = null;
  private byte[] mSdfQualities = null;
  private byte[] mQualities = new byte[CgUtils.CG_RAW_READ_LENGTH];
  private String mXQField;
  private Integer mAlignmentScoreAttr = null;
  private int mAlignmentScore;
  private char mPreviousChar;

  private boolean mIsValid;
  private String mInvalidReason;

  private final int mUnknownNt = DNA.N.ordinal();

  private boolean mIsReverseCompliment;

  private final int mUnknownsPenalty;

  /**
   * @param unknownsPenalty score penalty for unknown nucleotides
   */
  public SuperCigarValidator(int unknownsPenalty) {
    mUnknownsPenalty = unknownsPenalty;
  }

  /**
   * Expand a SAM record's quality field using Super Cigar quality overlap field.
   * @param qualityBuffer a 35-long byte array buffer for the expanded quality
   * @param qualityOverlap the overlap quality values
   * @param samQualities the flattened qualities from the SAM record
   * @param first true if first of pair
   * @param reverseCompliment true if reverse compliment
   * @param phredifyQualityOverlap true if the qualities in the overlap region should be converted to phred scale (i.e. have 33 subtracted)
   * @return the expanded quality array
   * @throws BadSuperCigarException if the <code>samQualities</code> length plus overlap length does not equal the expected value
   */
  public static byte[] expandCgSuperCigarQualities(byte[] samQualities, byte[] qualityBuffer, String qualityOverlap, boolean first, boolean reverseCompliment, boolean phredifyQualityOverlap) throws BadSuperCigarException {
    assert qualityBuffer.length == CgUtils.CG_RAW_READ_LENGTH;
    final int overlapQualityOffset = phredifyQualityOverlap ? FastaUtils.PHRED_LOWER_LIMIT_CHAR : 0;
    final int xqlength = qualityOverlap == null ? 0 : qualityOverlap.length();
    if (samQualities.length + xqlength != qualityBuffer.length) {
      throw new BadSuperCigarException("SAM record qualities plus " + SamUtils.CG_OVERLAP_QUALITY + " not expected length. Was: " + (samQualities.length + xqlength) + " expected: " + qualityBuffer.length);
    }
    if ((first && !reverseCompliment) || (reverseCompliment && !first)) {
      System.arraycopy(samQualities, 0, qualityBuffer, 0, CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE);
      for (int i = 0; i < xqlength; i++) {
        qualityBuffer[CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE + i] = (byte) (qualityOverlap.charAt(i) - overlapQualityOffset);
      }
      System.arraycopy(samQualities, CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE, qualityBuffer, CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE + xqlength, samQualities.length - CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE);
    } else {
      System.arraycopy(samQualities, 0, qualityBuffer, 0, samQualities.length - CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE);
      for (int i = 0; i < xqlength; i++) {
        qualityBuffer[samQualities.length - CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE + i] = (byte) (qualityOverlap.charAt(i) - overlapQualityOffset);
      }
      System.arraycopy(samQualities, samQualities.length - CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE, qualityBuffer, qualityBuffer.length - CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE, CgGotohEditDistance.CG_FINAL_FRAGMENT_SIZE);
    }
    if (reverseCompliment) {
      Utils.reverseInPlace(qualityBuffer);
    }
    return qualityBuffer;
  }

  void setData(SAMRecord samRecord, byte[] sdfRead, byte[] sdfQualities) throws BadSuperCigarException {
    mSamRecord = samRecord;
    mSdfRead = sdfRead;
    mSdfQualities = sdfQualities;
    assert mSdfRead.length == CgGotohEditDistance.CG_RAW_READ_LENGTH;
    assert mSdfQualities.length == CgGotohEditDistance.CG_RAW_READ_LENGTH;

    super.setCigar(samRecord.getStringAttribute(SamUtils.CG_SUPER_CIGAR), samRecord.getStringAttribute(SamUtils.CG_READ_DELTA));
    setTemplateStart(samRecord.getAlignmentStart() - 1);
    mIsReverseCompliment = samRecord.getReadNegativeStrandFlag();
    mXQField = samRecord.getStringAttribute(SamUtils.CG_OVERLAP_QUALITY);

    mQualities = expandCgSuperCigarQualities(samRecord.getBaseQualities(), mQualities, mXQField, samRecord.getFirstOfPairFlag(), mIsReverseCompliment, true);

    if (samRecord.hasAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE)) {
      mAlignmentScoreAttr = samRecord.getIntegerAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE);
    } else {
      mAlignmentScoreAttr = null;
    }
  }

  private void reset() {
    mIsValid = true;
    mInvalidReason = null;
    mAlignmentScore = 0;
  }

  /**
   * @return true if this record is valid
   */
  public boolean isValid() {
    return mIsValid;
  }
  /**
   * @return the human readable reason for this record being invalid, or null if valid
   */
  public String getInvalidReason() {
    return mInvalidReason;
  }

  private void setInvalid(String reason) {
    mIsValid = false;
    mInvalidReason = reason;
  }

  /**
   * @throws BadSuperCigarException if the cigar is not valid
   */
  @Override
  public void parse() throws BadSuperCigarException {
    reset();
    super.parse();
    if (mIsValid) {
      checkQualities();
      if (mAlignmentScoreAttr != null) {
        checkAlignmentScore();
      }
      assert getReadLength() != 0 || mReadDeltaPos == mReadDelta.length : "readDelta.len=" + mReadDelta.length + " but should be " + mReadDeltaPos;
    }
  }

  @Override
  protected void doChunk(char ch, int count) throws BadSuperCigarException {
    if (mIsValid) {
      mPreviousChar = '`';
      super.doChunk(ch, count);
    }
  }

  @Override
  protected byte getReadDelta(int pos) {
    if (pos >= mReadDelta.length && getReadPosition() != mSdfRead.length) {
      setInvalid("Read delta (" + SamUtils.CG_READ_DELTA + ") too short, " + mSamRecord.getSAMString().trim());
      return 0;
    }
    return mReadDelta[pos];
  }

  private byte getSdfReadNt() {
    byte nt;
//    do {
      if ((mIsReverseCompliment && (mSdfRead.length - 1 - getReadPosition()) < 0)
          || !mIsReverseCompliment && getReadPosition() >= mSdfRead.length) {
        throw new IllegalStateException("Read too short for actions, " + mSamRecord.getSAMString().trim());
      }
      nt = mIsReverseCompliment ? mSdfRead[mSdfRead.length - 1 - getReadPosition()] : mSdfRead[getReadPosition()];
      /*if (nt == DnaUtils.CG_SPACER_VALUE) {
        mReadPos++;
      } else*/ if (mIsReverseCompliment) {
        nt = DNA.complement(nt);
      }
//    } while (nt == DnaUtils.CG_SPACER_VALUE);
    return nt;
  }

  @Override
  protected void doReadHardClip() {
  }

  @Override
  protected void doReadSoftClip(int readNt) {
    if (!mIsValid) {
      return;
    }
    final int originalReadNt = getSdfReadNt();
    if (readNt != originalReadNt) {
      setInvalid("SDF read soft clip: " + DnaUtils.getBase(originalReadNt) + " does not match SAM " + SamUtils.CG_READ_DELTA + ": " + DnaUtils.getBase(readNt) + ", " + mSamRecord.getSAMString().trim());
    }
    mAlignmentScore++;
  }

  @Override
  protected void doTemplateOverlap() {
    if (!mIsValid) {
      return;
    }
//    System.err.println("Template overlap...");
    if (mXQField == null) {
      setInvalid("Overlap described but no " + SamUtils.CG_OVERLAP_QUALITY + " field present in SAM record, " + mSamRecord.getSAMString().trim());
    }
  }

  @Override
  protected void doReadOnly(int readNt) {
    if (!mIsValid) {
      return;
    }
    final int originalReadNt = getSdfReadNt();
    if (readNt != originalReadNt) {
      setInvalid("SDF read insert: " + DnaUtils.getBase(originalReadNt) + " does not match SAM " + SamUtils.CG_READ_DELTA + ": " + DnaUtils.getBase(readNt) + ", " + mSamRecord.getSAMString().trim());
    }
    mAlignmentScore += (mPreviousChar == 'I') ? 1 : 2;
    mPreviousChar = 'I';
  }

  @Override
  protected void doTemplateOnly(int templateNt) {
    mAlignmentScore += (mPreviousChar == 'D') ? 1 : 2;
    mPreviousChar = 'D';
  }

  @Override
  protected void doSubstitution(int readNt, int templateNt) {
    if (!mIsValid) {
      return;
    }
    final int originalReadNt = getSdfReadNt();

    if (originalReadNt == templateNt) {
      setInvalid("Expected mismatch, " + mSamRecord.getSAMString().trim());
    } else if (readNt != originalReadNt) {
      setInvalid(SamUtils.CG_READ_DELTA + " value: " + DnaUtils.getBase(readNt) + " does not match read value: " + DnaUtils.getBase(originalReadNt) + ", " + mSamRecord.getSAMString().trim());
    }
    mAlignmentScore++;
  }

  @Override
  protected void doEquality(int readNt, int nt) {
    if (!mIsValid) {
      return;
    }
    final int originalReadNt = getSdfReadNt();
    if (nt != mUnknownNt) {
      if (nt != originalReadNt && originalReadNt != mUnknownNt) {
        setInvalid("Expected match, SDF read=" + DnaUtils.getBase(originalReadNt) + ", template=" + DnaUtils.getBase(nt) + ", " + mSamRecord.getSAMString().trim());
      }
    }
  }

  @Override
  protected void doUnknownOnTemplate(int readNt, int templateNt) {
    if (!mIsValid) {
      return;
    }
    if (templateNt != mUnknownNt) {
      setInvalid("Expected unknown on template, was " + DnaUtils.getBase(templateNt));
    }
    final int originalReadNt = getSdfReadNt();
    if (readNt != originalReadNt) {
      setInvalid(SamUtils.CG_READ_DELTA + " value: " + DnaUtils.getBase(readNt) + " does not match read value: " + DnaUtils.getBase(originalReadNt) + ", " + mSamRecord.getSAMString().trim());
    }
    mAlignmentScore += mUnknownsPenalty;
  }

  @Override
  protected void doUnknownOnRead() {
    if (!mIsValid) {
      return;
    }
    if (getSdfReadNt() != mUnknownNt) {
      setInvalid("Expected unknown on read, was " + DnaUtils.getBase(getSdfReadNt()));
    }
    mAlignmentScore += mUnknownsPenalty;
  }

  private void checkQualities() {
    if (!Arrays.equals(mSdfQualities, mQualities)) {
      setInvalid("SDF and SAM qualities don't match, " + mSamRecord.getSAMString().trim());
    }

    /*for (int i = 0; i < mSdfQualities.length; i++) {
      if (mSdfQualities[i] != mQualities[i] + 33) {
        setInvalid("SDF and Sam qualities don't match, " + mSamRecord.format());
        return;
      }
    }*/
  }

  private void checkAlignmentScore() {
    if (mAlignmentScore != mAlignmentScoreAttr) {
      setInvalid("Super cigar alignment score was " + mAlignmentScore + ", but AS attribute was " + mAlignmentScoreAttr);
    }
  }
}
