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
package com.rtg.alignment;

import static com.rtg.alignment.ActionsHelper.CG_GAP_IN_READ;
import static com.rtg.alignment.ActionsHelper.CG_OVERLAP_IN_READ;
import static com.rtg.alignment.ActionsHelper.DELETION_FROM_REFERENCE;
import static com.rtg.alignment.ActionsHelper.INSERTION_INTO_REFERENCE;
import static com.rtg.alignment.ActionsHelper.MISMATCH;
import static com.rtg.alignment.ActionsHelper.NOOP;
import static com.rtg.alignment.ActionsHelper.SAME;
import static com.rtg.alignment.ActionsHelper.SOFT_CLIP;
import static com.rtg.alignment.ActionsHelper.UNKNOWN_READ;
import static com.rtg.alignment.ActionsHelper.UNKNOWN_TEMPLATE;

import com.reeltwo.jumble.annotations.JumbleIgnore;
import com.rtg.mode.DnaUtils;
import com.rtg.ngs.tempstage.BinaryTempFileRecord;
import com.rtg.pairedend.InsertHelper;
import com.rtg.sam.CigarFormatter;
import com.rtg.sam.MismatchPositionFormatter;
import com.rtg.sam.SamBamConstants;
import com.rtg.sam.SuperCigar;
import com.rtg.util.Utils;

/**
 * Store an alignment result and convert it to a SAM record.
 *
 */
public class AlignmentResult {

  private final byte[] mRead;
  private final byte[] mTemplate;
  private int mReadId;
  int mReferenceId;
  boolean mFirst;
  boolean mIsReverse = false;
  private final int mScore;
  private final int[] mActions;

  /**
   * True if <code>mActions</code> contains CG-specific commands.
   * This is true when the Probabilistic CG aligner was used.
   */
  private final boolean mCgGotohAligned;

  /**
   * Construct an Alignment Result object.
   * @param read bytes of the read
   * @param res input array of packed alignment actions
   * @param template bytes of the template
   */
  public AlignmentResult(byte[] read, int[] res, byte[] template) {
    mRead = read;
    mTemplate = template;

    assert res != null : "Null actions array passed into alignment result";
    mActions = ActionsHelper.copy(res);
    mCgGotohAligned = ActionsHelper.isCg(mActions);
    mScore = ActionsHelper.alignmentScore(mActions);
  }

  String cgReadString(StringBuilder sbRead) {
    //TODO: cg doesn't support SOFTCLIP or NOOP actions yet (they never get generated for cg)
    final ActionsHelper.CommandIterator iter;
    final int direction;
    int readPos;
    if (mFirst) {
      iter = ActionsHelper.iterator(mActions);
      direction = +1;
      readPos = 0;
    } else {
      iter = ActionsHelper.iteratorReverse(mActions);
      direction = -1;
      readPos = mRead.length - 1;
    }
    int cgOverlapSize = 0;
    while (iter.hasNext()) {
      final int action = iter.next();
      if (cgOverlapSize > 0 && action != CG_OVERLAP_IN_READ) {
        // start/continue skipping over cgOverlapSize nucleotides of the TEMPLATE after the overlap.
        if (action != INSERTION_INTO_REFERENCE) {
          --cgOverlapSize; // count down in template positions
        }
        if (action != DELETION_FROM_REFERENCE && action != CG_GAP_IN_READ) {
          readPos += direction;
        }
        continue;
      }
      switch (action) {
        case SAME:
        case INSERTION_INTO_REFERENCE:
        case MISMATCH:
        case UNKNOWN_TEMPLATE:
        case UNKNOWN_READ:
          sbRead.append(DnaUtils.getBase(mRead[readPos]));
          readPos += direction;
          break;
        case CG_OVERLAP_IN_READ:
          ++cgOverlapSize;
          break;
        case CG_GAP_IN_READ:
        case DELETION_FROM_REFERENCE:
          break;
        default:
          throw new RuntimeException("Illegal actions array: " + ActionsHelper.toString(mActions));
      }
    }
    assert readPos == (direction == 1 ? mRead.length : -1) : ActionsHelper.toString(mActions) + " but readPos=" + readPos;
    return mFirst ? sbRead.toString() : sbRead.reverse().toString();
  }

  /** @return the sanitised read string, with CG overlap bases and gap regions removed. */
  String readString() {
    final StringBuilder sbRead = new StringBuilder(mRead.length + 10);
    if (mCgGotohAligned) {
      return cgReadString(sbRead);
    } else {
      final ActionsHelper.CommandIterator iter = ActionsHelper.iterator(mActions);
      int r = 0;
      while (iter.hasNext()) {
        switch (iter.next()) {
          case DELETION_FROM_REFERENCE:
          case NOOP:
            break;
          default: // SAME, MISMATCH, INSERT, SOFTCLIP
            sbRead.append(DnaUtils.getBase(mRead[r++]));
            break;
        }
      }
    }
    return sbRead.toString();
  }

  /**
   * Preset the stats to keep the first/reverse info with the record, instead of
   * trying to infer it later
   * @param first whether is first in read pair
   * @param reverse whether is mapped on reverse strand
   */
  public void setIdentifyingInfo(boolean first, boolean reverse) {
    mFirst = first;
    mIsReverse = reverse;
  }

  /**
   * Set the rest of the output parameters if others already set
   * @param readId id of read
   * @param templateId name of template
   */
  public void setRemainingOutput(int readId, int templateId) {
    mReferenceId = templateId;
    mReadId = readId;
  }

  /** @return true if the hit is for the first read arm. */
  public boolean isFirst() {
    return mFirst;
  }

  /** @return the zero-based start position of the read after alignment */
  public int getStart() {
    return mActions[ActionsHelper.TEMPLATE_START_INDEX];
  }

  /** @return the read id for this alignment result */
  public int getReadId() {
    return mReadId;
  }

  /** @return the actions string */
  public String getActionsString() {
    return ActionsHelper.toString(mActions);
  }

  /** @return the number of matches in this alignment */
  public int getMatchCount() {
    return ActionsHelper.matchCount(mActions);
  }

  /** @return the number of deletions in the read (aka insertions into template) of this alignment */
  public int getDeletionsFromReadCount() {
    return ActionsHelper.deletionFromReadAndOverlapCount(mActions);
  }

  /** @return true if reverse */
  public boolean isReverse() {
    return mIsReverse;
  }

  /** @return the template name for this alignment result */
  public int getReferenceId() {
    return mReferenceId;
  }

  /** @return a simple count of the number of mismatches within aligned bases (i.e. excluding off-template or clipped) */
  protected int mismatches() {
    final ActionsHelper.CommandIterator iter = mIsReverse ? ActionsHelper.iteratorReverse(mActions) : ActionsHelper.iterator(mActions);
    int mismatches = 0;
    int tempPos = getStart();

    while (iter.hasNext() && tempPos < 0) {
      final int action = iter.next();
      switch (action) {
        case SAME:
        case MISMATCH:
        case UNKNOWN_TEMPLATE:
        case UNKNOWN_READ:
        case DELETION_FROM_REFERENCE:
        case CG_GAP_IN_READ:
          ++tempPos;
          break;
        case INSERTION_INTO_REFERENCE:
        case SOFT_CLIP:
          break;
        case CG_OVERLAP_IN_READ:
          --tempPos;
          break;
        default:
      }
    }
    while (iter.hasNext() && tempPos < mTemplate.length) {
      final int action = iter.next();
      if (action == DELETION_FROM_REFERENCE || action == MISMATCH || action == INSERTION_INTO_REFERENCE) {
        ++mismatches;
      }
      switch (action) {
        case SAME:
        case MISMATCH:
        case UNKNOWN_TEMPLATE:
        case UNKNOWN_READ:
        case DELETION_FROM_REFERENCE:
        case CG_GAP_IN_READ:
          ++tempPos;
          break;
        case INSERTION_INTO_REFERENCE:
        case SOFT_CLIP:
          break;
        case CG_OVERLAP_IN_READ:
          --tempPos;
          break;
        default:
      }
    }

    return mismatches;
  }

  String getCigarString(boolean rc, boolean legacy) {
    return CigarFormatter.actionsToCigar(mActions, rc, mTemplate.length, legacy, mFirst);
  }

  private void toRecord(BinaryTempFileRecord rec, boolean pairedEnd, int templateOffset, boolean legacy) {
    rec.setReadId(mReadId);
    rec.setReferenceId(mReferenceId);
    int samFlags = 0;
    if (mIsReverse) {
      samFlags |= SamBamConstants.SAM_READ_IS_REVERSE;
    }
    if (pairedEnd) {
      samFlags |= SamBamConstants.SAM_READ_IS_PAIRED;
      if (mFirst) {
        samFlags |= SamBamConstants.SAM_READ_IS_FIRST_IN_PAIR;
      } else {
        samFlags |= SamBamConstants.SAM_READ_IS_SECOND_IN_PAIR;
      }
    }
    rec.setSamFlags((byte) samFlags);
    // Reads matching off the left end of the template have CIGARs which are
    // soft-clipped on the left to make them start precisely at base 1 of the template.
    rec.setStartPosition(Math.max(1, ActionsHelper.zeroBasedTemplateStart(mActions) + templateOffset + 1));
    rec.setAlignmentScore(getScore());
    if (!mCgGotohAligned) {
      rec.setNumberMismatches(mismatches());
    }
    final String cigarString = getCigarString(mIsReverse, legacy);
    rec.setCigarString(cigarString.getBytes());
    if (legacy) {
      final String mismatchString = mCgGotohAligned ? "" : MismatchPositionFormatter.actionsToMismatchPositions(mActions, mIsReverse, mTemplate);
      if (mismatchString.length() > 0) {
        rec.setMdString(mismatchString.getBytes());
      } else {
        rec.setMdString(new byte[0]);
      }
    } else {
      rec.setMdString(new byte[0]);
    }
  }

  private void addCgGotohInfo(BinaryTempFileRecord rec) {
    rec.setSuperCigarString(SuperCigar.actionsToSuperCigar(mActions, mIsReverse, mTemplate.length).getBytes());
    final String delta = SuperCigar.readDelta(mActions, mRead, mTemplate.length);
    rec.setReadDeltaString(mIsReverse ? DnaUtils.reverseComplement(delta).getBytes() : delta.getBytes());
  }

  private void addMateInfo(BinaryTempFileRecord rec, AlignmentResult mateResult, int templateOffset) {
    if (mateResult != null) {
      final int thisStart = ActionsHelper.zeroBasedTemplateStart(mActions) + templateOffset;
      final int thisLength = ActionsHelper.templateLength(mActions);
      final int mateStart = ActionsHelper.zeroBasedTemplateStart(mateResult.mActions) + templateOffset;
      final int mateLength = ActionsHelper.templateLength(mateResult.mActions);
      final int signedTemplateLength = InsertHelper.tlen(mFirst, thisStart, thisLength, mateStart, mateLength);
      addMateInfo(rec, ActionsHelper.zeroBasedTemplateStart(mateResult.mActions), mateResult.mIsReverse, signedTemplateLength, templateOffset);
    }
  }

  private void addMateInfo(BinaryTempFileRecord rec, int zeroBasedMateStart, boolean materc,
      int inferredInsertSize,
      int templateOffset) {
    rec.setMatePosition(Math.max(1, zeroBasedMateStart + templateOffset + 1));
    rec.setTemplateLength(inferredInsertSize);
    int samFlags = rec.getSamFlags();
    samFlags |= SamBamConstants.SAM_READ_IS_MAPPED_IN_PROPER_PAIR;
    if (materc) {
      samFlags |= SamBamConstants.SAM_MATE_IS_REVERSE;
    }
    rec.setSamFlags((byte) samFlags);
    rec.setReferenceId(mReferenceId);
  }

  private void addCgAndQualities(BinaryTempFileRecord rec, String readString) {
    rec.setCgReadString(readString.getBytes());
    if (mCgGotohAligned) { //always add super cigar for anything aligned by CG Gotoh
      addCgGotohInfo(rec);
    }
  }

  /**
   * Convert the alignment to a record.
   *
   * @param pairedEnd are we in paired end mode
   * @param mateResult the AlignmentResult of the pair mate (null if none or single ended reads)
   * @param templateOffset template offset position
   * @param unfiltered true if we are running in <code>--allhits</code> (unfiltered) mode
   * @param legacy true iff using legacy cigars.
   * @return SAM record
   */
  public BinaryTempFileRecord toRecord(boolean pairedEnd, AlignmentResult mateResult, int templateOffset, boolean unfiltered, boolean legacy) {
    final BinaryTempFileRecord ret = new BinaryTempFileRecord(pairedEnd, legacy, mCgGotohAligned, unfiltered);

    toRecord(ret, pairedEnd, templateOffset, legacy);
    if (pairedEnd) {
      addMateInfo(ret, mateResult, templateOffset);
    }
    if (mCgGotohAligned) {
      final String readString = mIsReverse ? DnaUtils.reverseComplement(readString()) : readString();
      addCgAndQualities(ret, readString);
    }
    return ret;
  }

  /** @return the alignment score */
  public int getScore() {
    return mScore;
  }

  /**
   * Checks the template name and start position to see if this is the same result
   * Note: it is assumed that both have the same read id, and both are of this class
   * @param o object to compare against
   * @return whether considered equal
   */
  @Override
  public boolean equals(final Object o) {
    if (o == null) {
      return false;
    }
    final AlignmentResult ar = (AlignmentResult) o;
    return getStart() == ar.getStart() && mIsReverse == ar.mIsReverse && mReferenceId == ar.mReferenceId;
  }

  @Override
  @JumbleIgnore
  public int hashCode() {
    return Utils.pairHash(Utils.pairHash(getStart(), mIsReverse ? "R".hashCode() : "F".hashCode()), mReferenceId);
  }

  @Override
  public String toString() {
    return mReferenceId + "\t" + (mIsReverse ? "R" : "F") + "\t" + mReadId + "\t" + getStart() + "\t" + mScore + "\t" + ActionsHelper.toString(mActions) + ((!mCgGotohAligned) ? "\t" + mismatches() : "");
  }
}
