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
  private byte[] mQualities = null;
  private String mXQField;
  private Integer mAlignmentScoreAttr = null;
  private int mAlignmentScore;
  private char mPreviousChar;

  private boolean mIsValid;
  private String mInvalidReason;

  private final int mUnknownNt = DNA.N.ordinal();

  private boolean mIsReverseComplement;

  private final int mUnknownsPenalty;

  //keeps track of whether we are currently replaying a template position
  //-ve means we are reusing template
  //+ve means we are skipping template
  //The basic idea is to check whether an overlap has at least as many deletions
  //immediately surrounding it.
  //this doesn't account for certain nasty theoretical alignments however.
  private int mOverlapTemplate;

  private boolean mReusedTemplate;


  /**
   * @param unknownsPenalty score penalty for unknown nucleotides
   */
  public SuperCigarValidator(int unknownsPenalty) {
    mUnknownsPenalty = unknownsPenalty;
  }

  /**
   * Expand a SAM record's quality field using Super Cigar quality overlap field. The result will still be
   * in SAM orientation, so may require reversing if the original read orientation is required.
   * @param dest the destination buffer for the expanded quality, which must be the correct length for the CG read type.
   * @param samQualities the flattened qualities from the SAM record
   * @param qualityOverlap the overlap quality values (<code>XQ</code> field)
   * @param first true if the read is first of pair
   * @param reverse true if the alignment was in the the reverse orientation
   * @throws BadSuperCigarException if the <code>samQualities</code> length plus overlap length does not equal the expected value
   */
  public static void unrollSuperCigarQualities(byte[] dest, byte[] samQualities, String qualityOverlap, boolean first, boolean reverse) throws BadSuperCigarException {
    final int xqlength = qualityOverlap == null ? 0 : qualityOverlap.length();
    final byte[] xq = qualityOverlap == null ? null : FastaUtils.asciiToRawQuality(qualityOverlap);
    final int unrollLen = samQualities.length + xqlength;
    if (unrollLen != CgUtils.CG_RAW_READ_LENGTH && unrollLen != CgUtils.CG2_RAW_READ_LENGTH) {
      throw new BadSuperCigarException("SAM record qualities plus " + SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY + " not expected length. Was: " + unrollLen + " expected: " + CgUtils.CG_RAW_READ_LENGTH + " or " + CgUtils.CG2_RAW_READ_LENGTH);
    }
    assert dest.length == CgUtils.CG_RAW_READ_LENGTH || dest.length == CgUtils.CG2_RAW_READ_LENGTH;
    final boolean v2 = dest.length == CgUtils.CG2_RAW_READ_LENGTH;
    final int fwdOverlap = v2 ? CgUtils.CG2_OVERLAP_POSITION : CgUtils.CG_OVERLAP_POSITION;
    final boolean overlapOnLeft = v2 && !reverse || !v2 && first == !reverse;
    final int overlapPos = overlapOnLeft ? fwdOverlap : samQualities.length - fwdOverlap;
    System.arraycopy(samQualities, 0, dest, 0, overlapPos);
    if (xq != null) {
      System.arraycopy(xq, 0, dest, overlapPos, xq.length);
    }
    System.arraycopy(samQualities, overlapPos, dest, overlapPos + xqlength, samQualities.length - overlapPos);
  }

  void setData(SAMRecord samRecord, byte[] sdfRead, byte[] sdfQualities) throws BadSuperCigarException {
    mSamRecord = samRecord;
    mSdfRead = sdfRead;
    mSdfQualities = sdfQualities;
    assert mSdfRead.length == CgUtils.CG_RAW_READ_LENGTH || mSdfRead.length == CgUtils.CG2_RAW_READ_LENGTH;
    assert mSdfQualities.length == mSdfRead.length;

    super.setCigar(samRecord.getStringAttribute(SamUtils.CG_SUPER_CIGAR), samRecord.getStringAttribute(SamUtils.CG_READ_DELTA));
    setTemplateStart(samRecord.getAlignmentStart() - 1);
    mIsReverseComplement = samRecord.getReadNegativeStrandFlag();
    mXQField = samRecord.getStringAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY);

    if (mQualities == null || mQualities.length != mSdfQualities.length) {
      mQualities = new byte[mSdfQualities.length];
    }
    unrollSuperCigarQualities(mQualities, samRecord.getBaseQualities(), mXQField, samRecord.getFirstOfPairFlag(), mIsReverseComplement);
    if (mIsReverseComplement) {
      Utils.reverseInPlace(mQualities);
    }

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
    mOverlapTemplate = 0;
    mReusedTemplate = false;
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
      if (mIsValid) {
        checkOverlap();
      }
      assert getReadLength() != 0 || mReadDeltaPos == mReadDelta.length : "readDelta.len=" + mReadDelta.length + " but should be " + mReadDeltaPos;
    }
  }

  protected void checkOverlap() {
    if (mReusedTemplate && mXQField == null) {
      setInvalid("Overlap described but no " + SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY + " field present in SAM record, " + mSamRecord.getSAMString().trim());
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
      if ((mIsReverseComplement && (mSdfRead.length - 1 - getReadPosition()) < 0)
          || !mIsReverseComplement && getReadPosition() >= mSdfRead.length) {
        throw new IllegalStateException("Read too short for actions, " + mSamRecord.getSAMString().trim());
      }
      nt = mIsReverseComplement ? mSdfRead[mSdfRead.length - 1 - getReadPosition()] : mSdfRead[getReadPosition()];
      /*if (nt == DnaUtils.CG_SPACER_VALUE) {
        mReadPos++;
      } else*/ if (mIsReverseComplement) {
        nt = DNA.complement(nt);
      }
//    } while (nt == DnaUtils.CG_SPACER_VALUE);
    return nt;
  }

  @Override
  protected void doReadHardClip() {
    templateOnlyOp();
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

    mOverlapTemplate--;
//    System.err.println("Template overlap...");
//    if (mXQField == null) {
//      setInvalid("Overlap described but no " + SamUtils.CG_OVERLAP_QUALITY + " field present in SAM record, " + mSamRecord.getSAMString().trim());
//    }
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
    templateOnlyOp();
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
    readAndTemplate();
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
    readAndTemplate();
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

  @Override
  protected void doUnmapped() {
    super.doUnmapped();
    readAndTemplate();
  }

  @Override
  protected void doTemplateSkip(int templateNt) {
    super.doTemplateSkip(templateNt);
    templateOnlyOp();
  }

  private void readAndTemplate() {
    if (mOverlapTemplate < 0) {
      mReusedTemplate = true;
      mOverlapTemplate++;
    }
  }

  private void templateOnlyOp() {
    mOverlapTemplate++;
  }
}
