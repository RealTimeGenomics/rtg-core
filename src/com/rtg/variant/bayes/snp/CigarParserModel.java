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
package com.rtg.variant.bayes.snp;


import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.mode.DNA;
import com.rtg.reader.CgUtils;
import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.SamUtils;
import com.rtg.sam.SuperCigarParser;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.ReadParserInterface;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;
import com.rtg.variant.util.VariantUtils;

/**
 * Generate matcher reports from a CIGAR.
 */
public class CigarParserModel implements ReadParserInterface {

  private static final boolean SOFT_CLIP_COMPLEX_TRIGGER = GlobalFlags.getBooleanValue(CoreGlobalFlags.SOFT_CLIP_COMPLEX_TRIGGER);

  /** Handles explicit nucleotides and deletions. */
  private final MatcherInterface mMatcherNt;

  /** Handles inserts. */
  private final MatcherInterface mMatcherIndel;

  private final int mStart;

  private final int mEnd;

  private final ToMatcherCigarParser mParser = new ToMatcherCigarParser();

  private final VariantParams mParams;

  /**
   * Construct a new formatter using the specified matchers.
   * @param matcherNt matcher for nucleotides and deletes.
   * @param indelMatcher matcher for insertions (can be null).
   * @param start of region to be processed (0 based, inclusive)
   * @param end of region to be processed (0 based, exclusive)
   * @param params command line parameters.
   */
  public CigarParserModel(final MatcherInterface matcherNt, MatcherInterface indelMatcher, int start, int end, VariantParams params) {
    mMatcherNt = matcherNt;
    mMatcherIndel = indelMatcher;
    mStart = start;
    mEnd = end;
    mParams = params;
  }

  int getReadScore(VariantAlignmentRecord var) {
    return VariantUtils.readScoreFromAlignmentRecord(var, mParams);
  }

  @Override
  public void toMatcher(VariantAlignmentRecord var, MachineType machineType, int qdefault, byte[] templateBytes) throws BadSuperCigarException {
    final byte[] qualities;
    try {
      // Note that we don't use the SuperCigar for CG here, as the overlaps would be treated as double-evidence (whereas allpaths does the right thing).
      mParser.setStandardCigar(var.getCigar(), var.getRead(), var.getRead().length);
    } catch (final IllegalArgumentException iae) {
      throw new BadSuperCigarException("Illegal DNA character", iae);
    }
    qualities = var.getRecalibratedQuality().length == 0 || mParams.ignoreQualityScores() ? null : var.getRecalibratedQuality();

    mParser.setTemplateStart(var.getStart());
    mParser.setTemplate(templateBytes);
    mParser.setQualities(qualities);
    mParser.setAdditional(machineType, getReadScore(var), qdefault, var.isFirst(), !var.isNegativeStrand(), var.isReadPaired(), var.isMated(), var.isUnmapped());

    mParser.parse();
  }

  class ToMatcherCigarParser extends SuperCigarParser {

    ToMatcherCigarParser() { }

    private boolean mIsForward;
    private boolean mIsReadPaired;
    private boolean mIsMated;
    private byte[] mQualities = null;
    private boolean mCgTrimOuterBases;
    private boolean mIsCgOverlapLeft = false;
    private int mMapScore;
    private int mDefaultQuality;
    private boolean mInBlock;

    @Override
    public void setCigar(String superCigar, String readDelta) {
      super.setCigar(superCigar, readDelta);
      mQualities = null;
    }

    @Override
    public void setStandardCigar(String cigar, byte[] read, int readLength) {
      super.setStandardCigar(cigar, read, readLength);
      mQualities = null;
    }

    /**
     * Set the expanded qualities for this SAM record
     * @param qualities the expanded qualities
     */
    void setQualities(byte[] qualities) {
      mQualities = qualities;
    }

    int getCurrentQuality() throws BadSuperCigarException {
      if (mQualities != null && getReadPosition() >= mQualities.length) {
        throw new BadSuperCigarException("readPos " + getReadPosition() + " > qual.len in SAM record");
      }
      return mQualities == null ? mDefaultQuality : mQualities[getReadPosition()];
    }

    void setAdditional(MachineType machineType, int mapScore, int defaultQuality, boolean isFirst, boolean isForward, boolean isReadPaired, boolean isMated, boolean isUnmapped) {
      mCgTrimOuterBases = machineType != null && machineType.cgTrimOuterBases();
      mMapScore = mapScore;
      mDefaultQuality = defaultQuality;
      mIsCgOverlapLeft = machineType == MachineType.COMPLETE_GENOMICS && (isFirst ^ !isForward)
        || machineType == MachineType.COMPLETE_GENOMICS_2 && isForward;
      mIsForward = isForward;
      mIsReadPaired = isReadPaired;
      mIsMated = isMated;
      mIsUnmapped = isUnmapped;
    }

    @Override
    protected void doChunk(char ch, int count) throws BadSuperCigarException {
      mInBlock = false;
      if ((ch == SamUtils.CIGAR_SAME || ch == SamUtils.CIGAR_MISMATCH || ch == SamUtils.CIGAR_SAME_OR_MISMATCH || ch == SamUtils.CIGAR_DELETION_FROM_REF)
          && getTemplatePosition() + count > getTemplateLength()) {
        throw new BadSuperCigarException("Cigar exceeds template length");
      }
      super.doChunk(ch, count);
    }

    @Override
    public void parse() throws BadSuperCigarException {
      super.parse();
    }

    @Override
    protected void doReadOnly(int readNt) {
      if (mInBlock) {
        return;
      }
      mInBlock = true;
      final int templatePosition = getTemplatePosition();
      if (mMatcherIndel != null && includeBase(getReadPosition()) && mStart <= templatePosition && templatePosition < mEnd) {
        mMatcherIndel.match(templatePosition, EvidenceIndelFactory.SINGLETON.evidence(EvidenceIndel.INSERT, /*unused*/0, /*unused*/0, mMapScore, /*unused*/0, /*unused*/0, mOperationLength, false));
      }
    }

    @Override
    protected void doReadSoftClip(int readNt) {
      if (SOFT_CLIP_COMPLEX_TRIGGER) {
        if (mInBlock) {
          return;
        }
        mInBlock = true;
        final int readPosition = getReadPosition();
        final int templatePosition = getTemplatePosition();
        if (mMatcherIndel != null && includeBase(readPosition) && mStart <= templatePosition && templatePosition < mEnd) {
          mMatcherIndel.match(templatePosition, EvidenceIndelFactory.SINGLETON.evidence(readPosition == 0 ? EvidenceIndel.SOFT_CLIP_LEFT : EvidenceIndel.SOFT_CLIP_RIGHT, /*unused*/0, /*unused*/0, mMapScore, /*unused*/0, /*unused*/0, mOperationLength, false));
        }
      }
    }

    @Override
    protected void doTemplateOnly(int templateNt) {
      if (mMatcherIndel != null && includeBase(getReadPosition())) {
        final int templatePosition = getTemplatePosition();
        if (mStart <= templatePosition && templatePosition < mEnd) {
          mMatcherIndel.match(templatePosition, EvidenceIndelFactory.SINGLETON.evidence(EvidenceIndel.DELETE, /*unused*/0, /*unused*/0, mMapScore, /*unused*/0, /*unused*/0, mOperationLength, false));
        }
      }
    }

    @Override
    protected void doSubstitution(int readNt, int templateNt) throws BadSuperCigarException {
      assert readNt != templateNt : mCigar + "@rpos=" + getReadPosition() + " @tstartpos=" + getTemplateStartPosition() + " : r" + readNt + " == t" + templateNt;
      doSubstitutionOrEquality(readNt);
    }

    @Override
    protected void doEquality(int readNt, int nt) throws BadSuperCigarException {
      assert readNt == nt || nt == DNA.N.ordinal() || readNt == DNA.N.ordinal() : mCigar + "@pos=" + getReadPosition() + " @tstartpos=" + getTemplateStartPosition() + " : r" + readNt + " != t" + nt;
      doSubstitutionOrEquality(readNt);
    }

    private void doSubstitutionOrEquality(final int readNt) throws BadSuperCigarException {
      final int readPosition = getReadPosition();
      final boolean should = includeBase(readPosition);
      final int templatePosition = getTemplatePosition();
      if (should && mStart <= templatePosition && templatePosition < mEnd) {
        mMatcherNt.match(templatePosition, readPosition, getReadLength() - readPosition - 1, readNt, mMapScore, getCurrentQuality(), mMatcherNt.getStateIndex(mIsForward, mIsReadPaired, mIsMated));
        if (mMatcherIndel != null) { //for coverage-sensitive indel triggering, we create a dummy match for the indels to compute coverage down deeper.
          mMatcherIndel.match(templatePosition, null);
        }
      }
    }

    @Override
    protected void doUnmapped() {
      final int templatePosition = getTemplatePosition();
      if (mStart <= templatePosition && templatePosition < mEnd) {
        mMatcherNt.unmapped(templatePosition);
      }
    }

    @Override
    protected void doUnknownOnTemplate(int readNt, int templateNt) throws BadSuperCigarException {
      doSubstitutionOrEquality(readNt);
    }

    @Override
    protected void doUnknownOnRead() throws BadSuperCigarException {
      doSubstitutionOrEquality((byte) DNA.N.ordinal());
    }

    /**
     * Determine whether a base should be included, allowing the ability to trim the overlap end of CG version 1 reads
     * @param rPos read position
     * @return true if this base should be included
     */
    private boolean includeBase(int rPos) {
      final int readLength = getReadLength();
      return !mCgTrimOuterBases
        || (!mIsCgOverlapLeft && rPos < (CgUtils.CG_RAW_READ_LENGTH - CgUtils.CG_OVERLAP_POSITION))
        || (mIsCgOverlapLeft && rPos >= (readLength == 0 ? CgUtils.CG_OVERLAP_POSITION : readLength - (CgUtils.CG_RAW_READ_LENGTH - CgUtils.CG_OVERLAP_POSITION)));
    }

  }
}
