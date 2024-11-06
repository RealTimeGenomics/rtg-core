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
package com.rtg.variant.bayes.multisample;


import java.io.ByteArrayOutputStream;

import com.rtg.mode.DnaUtils;
import com.rtg.sam.SamUtils;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;
import com.rtg.variant.match.AlignmentMatch;
import com.rtg.variant.util.VariantUtils;

/**
 * Used to convert variant alignment records to alignment matches
 */
public class AlignmentRecordMatcher {

  private final VariantAlignmentRecord mVariantAlignmentRecord;
  private final String mCigar;
  private final byte[] mRead;
  private final byte[] mQual;
  private int mReadPos;
  private int mRefPos;

  private int mN;
  private int mI;
  private boolean mLeftN;
  private boolean mRightN;
  private boolean mValidNt;
  private boolean mReadNt;
  private int mSoftClippedStart;
  private int mSoftClippedEnd;
  private boolean mInvalid;

  private final StringBuilder mMatchString;
  private final ByteArrayOutputStream mQualString;

  /**
   * @param alignmentRecord alignment record to parse
   */
  public AlignmentRecordMatcher(VariantAlignmentRecord alignmentRecord) {
    mVariantAlignmentRecord = alignmentRecord;
    mCigar = alignmentRecord.getCigar();
    mRead = alignmentRecord.getRead();
    mQual = alignmentRecord.getRecalibratedQuality();
    //position in reference (0 based)
    mRefPos = alignmentRecord.getStart();
    //position in read (0 based)
    mReadPos = 0;
    mN = 0;
    mI = 0;
    mMatchString = new StringBuilder();
    mQualString = mQual.length > 0 ? new ByteArrayOutputStream() : null;
    mSoftClippedStart = 0;
    mSoftClippedEnd = mRead.length;
  }

  /**
   * Parses using match action, until reference position reaches given value or the end of cigar is reached
   * @param action function to invoke for each operation
   * @param untilPos stop parsing once reference reaches this position
   * @param greedy consume non reference matching operations at <code>untilPos</code>
   */
  private void parse(MatchAction action, int untilPos, boolean greedy) {
    while (mI < mCigar.length()) {
      final char c = mCigar.charAt(mI);
      if (Character.isDigit(c)) {
        mN = 10 * mN + c - '0';
      } else {
        while ((mN > 0) && (!greedy ? mRefPos < untilPos : greedyCheck(c, mRefPos, untilPos))) {
          action.match(c);
          postMatch(c);
          --mN;
        }
        if (mN > 0 && mRefPos >= untilPos) {
          break;
        }
      }
      ++mI;
    }
  }

  private static boolean greedyCheck(char op, int refPos, int untilPos) {
    return op == SamUtils.CIGAR_INSERTION_INTO_REF ? refPos <= untilPos : refPos < untilPos;
  }


  interface MatchAction {
    void match(char cigarOp);
  }

  final void postMatch(char c) {
    switch (c) {
      case SamUtils.CIGAR_SAME_OR_MISMATCH:
      case SamUtils.CIGAR_SAME:
      case SamUtils.CIGAR_MISMATCH:
        ++mRefPos;
        ++mReadPos;
        break;
      case SamUtils.CIGAR_INSERTION_INTO_REF:
      case SamUtils.CIGAR_SOFT_CLIP: // soft-clipping bases in read ignored for position
        ++mReadPos;
        break;
      case SamUtils.CIGAR_GAP_IN_READ:
      case SamUtils.CIGAR_DELETION_FROM_REF:
        ++mRefPos;
        break;
      case SamUtils.CIGAR_PADDING:
      case SamUtils.CIGAR_HARD_CLIP:
        break;
      default:
        throw new IllegalArgumentException(mCigar);
    }
  }

  final void startMatch(char c) {
    switch (c) {
      case SamUtils.CIGAR_SOFT_CLIP:
        adjustSoftClipBounds();
        break;
      case SamUtils.CIGAR_SAME_OR_MISMATCH:
      case SamUtils.CIGAR_SAME:
      case SamUtils.CIGAR_MISMATCH:
      case SamUtils.CIGAR_INSERTION_INTO_REF:
        mReadNt = true;
        break;
      case SamUtils.CIGAR_GAP_IN_READ:
      case SamUtils.CIGAR_DELETION_FROM_REF:
      case SamUtils.CIGAR_PADDING:
      case SamUtils.CIGAR_HARD_CLIP:
        break;
      default:
        throw new IllegalArgumentException(mCigar);
    }
  }

  final void middleMatch(char c) {
    switch (c) {
      case SamUtils.CIGAR_SAME_OR_MISMATCH:
      case SamUtils.CIGAR_SAME:
      case SamUtils.CIGAR_MISMATCH:
      case SamUtils.CIGAR_INSERTION_INTO_REF:
        if (mRightN || mReadPos >= mRead.length) {
          mInvalid = true;
        } else {
          mMatchString.append(DnaUtils.getBase(mRead[mReadPos]));
          if (mQualString != null) {
            mQualString.write(mQual[mReadPos]);
          }
          mValidNt = true;
        }
        break;
      case SamUtils.CIGAR_GAP_IN_READ:
        if (mValidNt) {
          mRightN = true;
        } else {
          mLeftN = true;
        }
        break;
      case SamUtils.CIGAR_SOFT_CLIP:
        adjustSoftClipBounds();
        break;
      case SamUtils.CIGAR_DELETION_FROM_REF:
      case SamUtils.CIGAR_PADDING:
      case SamUtils.CIGAR_HARD_CLIP:
        break;
      default:
        throw new IllegalArgumentException(mCigar);
    }
  }

  private void adjustSoftClipBounds() {
    if (mReadNt || mValidNt) {
      if (mReadPos < mSoftClippedEnd) {
        mSoftClippedEnd = mReadPos;
      }
    } else {
      mSoftClippedStart = mReadPos + 1;
    }
  }

  /**
   * Convert the alignment record into an alignment match by parsing its cigar. This method may only be called once.
   * @param chooser For scaling quality values
   * @param start start reference position (0-based inclusive)
   * @param end end reference position (0-based exclusive)
   * @param params variant calling parameters
   * @return the alignment match
   */
  public AlignmentMatch getMatch(MachineErrorChooserInterface chooser, int start, int end, VariantParams params) {
    if (mI != 0) {
      throw new RuntimeException("Not re-entrant");
    }
    if (mRead.length == 0) {
      return null; // For records that are non-primary and have no read or quality data stored with them
    }
    parse(this::startMatch, start, false);

    if (mI >= mCigar.length()) {
      return null;
    }
    final int startInRead = mReadPos;
    mLeftN = mRefPos > start;

    parse(this::middleMatch, end, true);

    if (mInvalid) {
      return null;
    }
    if (mI >= mCigar.length()) {
      mRightN = mRefPos < end;
    }
    final int endReadPos = mReadPos;
    //distinguish case when insert at start and start and end are the same from there being no insert at start
    if (mReadPos == 0 && mRefPos == start) {
      return null;
    }
    if (mLeftN && !mValidNt && !mRightN) {
      //all Ns - no useful match can be returned.
      return null;
    }
    final byte[] quality = mQualString == null ? null : mQualString.toByteArray();
    final String rs = mMatchString.toString();
    final AlignmentMatch match = new AlignmentMatch(mVariantAlignmentRecord, chooser, rs, quality, params.qDefault(), 0, rs.length(), VariantUtils.readScoreFromAlignmentRecord(mVariantAlignmentRecord, params), !mLeftN, !mRightN);
    match.setBasesLeftOfMatch(startInRead);
    match.setBasesRightOfMatch(mVariantAlignmentRecord.getRead().length - endReadPos);
    parse(this::startMatch, Integer.MAX_VALUE, false);
    match.setSoftClipLeft(mSoftClippedStart);
    match.setSoftClipRight(mSoftClippedEnd);

    return match;
  }
}
