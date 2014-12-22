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

package com.rtg.simulation;

import com.rtg.sam.SamUtils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Used by <code>ReadMappingAccuracy</code> to evaluate the start position
 * offset in the soft clipped region at the start of an alignment.
 */
public class SoftClipCigarParser {
  private final CigarActionState mGeneratedCigarState;
  private final CigarActionState mAlignedCigarState;

  /**
   * Constructs a parser which steps through a generated read's cigar
   * simultaneously with a cigar from an alignment, attempting to determine
   * if the start positions agree in the presence of soft clipping at the start
   * of the alignment.
   */
  public SoftClipCigarParser() {
    mGeneratedCigarState = new CigarActionState();
    mAlignedCigarState = new CigarActionState();
  }

  /**
   * Sets state of the parser. To be called per alignment checked.
   * @param generatedStartPos the zero based start position of the generated read on the template
   * @param alignedStartPos the zero based start position of the alignment on the template
   * @param generatedCigar the cigar of the generated read
   * @param alignedCigar the cigar of the alignment
   */
  void setState(long generatedStartPos, int alignedStartPos, String generatedCigar, String alignedCigar) {
    mGeneratedCigarState.init(generatedCigar, generatedStartPos);
    mAlignedCigarState.init(alignedCigar, alignedStartPos);
  }

  /**
   * Parse the two cigars.
   * @return the offset due to indels occurring in the region that was soft clipped
   */
  int parse() {
    int indelOffset = 0;
    char alignCigarChar = mAlignedCigarState.getAction();
    char genCigarChar = mGeneratedCigarState.getAction();

    if (alignCigarChar == SamUtils.CIGAR_SOFT_CLIP) {
      final int softClipLength = mAlignedCigarState.getCurrentActionLength();
      for (int i = 0; i < softClipLength; i++) {
        //count up the indels in the soft clipping region
        if (genCigarChar == SamUtils.CIGAR_DELETION_FROM_REF) {
          indelOffset--;
          i--; //we don't consume a soft clip for this case, hence the i--
        } else if (genCigarChar == SamUtils.CIGAR_INSERTION_INTO_REF) {
          indelOffset++;
          alignCigarChar = mAlignedCigarState.getAction();  //skip over soft clip
        } else {  //currently the only other type of action we generate is MNP or same.
          alignCigarChar = mAlignedCigarState.getAction();  //skip over soft clip
        }
        genCigarChar = mGeneratedCigarState.getAction();
      }
    } else {
      return 0; //no soft clip, no offset.
    }

    //we're out of the soft clip region.

    if (mGeneratedCigarState.getReadPos() != mAlignedCigarState.getReadPos()) { //we should be in sync readpos wise, if not, bail.
      return 0; // "noop" value
    } else if (mGeneratedCigarState.getTemplatePos() == mAlignedCigarState.getTemplatePos() && mGeneratedCigarState.getReadPos() == mAlignedCigarState.getReadPos()) {
      //we're actually in sync. we can stop
      return indelOffset;
    }

    //examine the next generated cigar action in case there's an indel right where the soft clipping occurs.
    if (genCigarChar == SamUtils.CIGAR_DELETION_FROM_REF) {
      indelOffset -= mGeneratedCigarState.getCurrentActionRemaining();
      if (alignCigarChar == SamUtils.CIGAR_DELETION_FROM_REF) { //bizarre case where the alignment is soft clip followed by a delete. No sane aligner would ever do this.
        //make it the difference between the two
        indelOffset += mAlignedCigarState.getCurrentActionRemaining();
      }
      return indelOffset; //look no further
    } else if (genCigarChar == SamUtils.CIGAR_INSERTION_INTO_REF) { //not actually sure if this case can occur
      indelOffset += mGeneratedCigarState.getCurrentActionRemaining();
      if (alignCigarChar == SamUtils.CIGAR_INSERTION_INTO_REF) { //bizarre case where the alignment is soft clip followed by an insert. No sane aligner would ever do this.
        //make it the difference between the two
        indelOffset -= mAlignedCigarState.getCurrentActionRemaining();
      }
      return indelOffset; //look no further
    }

    //now deal with the possible ambiguous indel which may or may not exist in generated.

    //start looking through the read bases at the two different positions. They should be the same, even though the template positions are different. If they aren't, this isn't an ambiguous region!
    while (genCigarChar != (char) -1 && alignCigarChar != (char) -1
        && (genCigarChar == alignCigarChar || alignCigarChar == SamUtils.CIGAR_SAME_OR_MISMATCH && (genCigarChar == SamUtils.CIGAR_SAME || genCigarChar == SamUtils.CIGAR_MISMATCH))) {
      //we've got two reads using the same bases at offsets on the template. The cigars should be the same, otherwise we're actually at different places.
      genCigarChar = mGeneratedCigarState.getAction();
      alignCigarChar = mAlignedCigarState.getAction();
    }
    if (genCigarChar == SamUtils.CIGAR_DELETION_FROM_REF) {
      indelOffset -= mGeneratedCigarState.getCurrentActionLength();
      //TODO for insurance, make sure the action and read/template positions are the same after this.
    } else if (genCigarChar == SamUtils.CIGAR_INSERTION_INTO_REF) { //not entirely sure this case is possible...
      indelOffset += mGeneratedCigarState.getCurrentActionLength();
      //TODO for insurance, make sure the action and read/template positions are the same after this.
    } //TODO doesn't REALLY make sense to see the indel in the aligned read here, but I suppose it MIGHT happen...

    return indelOffset;
  }

  static class CigarActionState {
    private char mCurrentAction;
    private int mCurrentActionLength;
    private int mCurrentActionRemaining;

    private String mCigar;
    private int mCigarPos;
    private int mReadPos;
    private long mTemplatePos;

    void init(String cigar, long templatePos) {
      mCigar = cigar;
      mReadPos = 0;
      mTemplatePos = templatePos;
      mCigarPos = 0;
      mCurrentActionLength = 0;
      mCurrentAction = (char) -1;
      mCurrentActionRemaining = 0;
    }

    /**
     * Get an action from the cigar. Note the read pos/template pos after this call will be for the NEXT base.
     * @return the next cigar action, or '-1' if we've reached the end of the cigar
     */
    char getAction() {
      if (mCurrentAction != (char) -1) {
        //deal with the previous action
        if (mCurrentAction == SamUtils.CIGAR_SAME || mCurrentAction == SamUtils.CIGAR_MISMATCH || mCurrentAction == SamUtils.CIGAR_SAME_OR_MISMATCH) {
          mReadPos++;
          mTemplatePos++;
        } else if (mCurrentAction == SamUtils.CIGAR_DELETION_FROM_REF) {
          mTemplatePos++;
        } else if (mCurrentAction == SamUtils.CIGAR_INSERTION_INTO_REF) {
          mReadPos++;
        } else if (mCurrentAction == SamUtils.CIGAR_SOFT_CLIP) {
          mReadPos++;
        } else if (mCurrentAction == SamUtils.CIGAR_GAP_IN_READ) {
          mTemplatePos++;
        } else {
          Diagnostic.developerLog("Unchecked action in getAction: " + mCurrentAction);
        }
        mCurrentActionRemaining--;
      }

      //get a new action if necessary
      if (mCurrentActionRemaining == 0 && !ensureAction()) {
        return (char) -1;
      }
      return mCurrentAction;
    }

    /**
     * Ensures an action is loaded from the cigar.
     * @return true if an action is available, or false if we've reached the end of the cigar
     */
    private boolean ensureAction() {
      if (mCigarPos == mCigar.length()) {
        return false;
      }
      if (mCigarPos == 0 || (mCigarPos < mCigar.length() && mCurrentActionRemaining == 0)) {  //haven't loaded anything, or are at end of previous action and there's more to go
        int mCigarCount = 0;
        while (mCigarPos < mCigar.length()) {
          final char ch = mCigar.charAt(mCigarPos);
          mCigarPos++;
          if (ch >= '0' && ch <= '9') {
            mCigarCount = mCigarCount * 10 + (ch - '0');
          } else {  //once we've found a non-numerical, break out
            mCurrentAction = ch;
            mCurrentActionLength = mCigarCount;
            mCurrentActionRemaining = mCigarCount;
            return true;
          }
        }
        return false;
      }
      return true;
    }

    int getReadPos() {
      return mReadPos;
    }

    long getTemplatePos() {
      return mTemplatePos;
    }

    int getCurrentActionLength() {
      return mCurrentActionLength;
    }

    int getCurrentActionRemaining() {
      return mCurrentActionRemaining;
    }
  }
}
