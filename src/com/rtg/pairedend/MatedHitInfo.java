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
package com.rtg.pairedend;


import com.rtg.alignment.AlignmentResult;

/**
 * A simple structure to store information about a pair of hits.
 * Left/right refers to placement on the reference.
 * First/second refers to whether the read arm was first or second in sequencing.
 *
 */
public class MatedHitInfo {

  int mReadId;
  boolean mFirstRight;
  boolean mReverseComplementRight;
  int mTemplateStartRight;
  boolean mReverseComplementLeft;
  int mTemplateStartLeft;

  // Populated and stored here so that alignments for a mating need only be calculated once
  AlignmentResult mAlignmentLeft = null;
  int mAlignmentScoreLeft = -1;
  AlignmentResult mAlignmentRight = null;
  int mAlignmentScoreRight = -1;

  // Used to indicate that the left alignment is so poor that there is no point attempting further matings with this hit
  boolean mLeftIsVeryPoor;

  /**
   * Helper method to set many values at once
   * @param readId the read id
   * @param firstRight true if the read that was first in sequencing is rightmost on the reference
   * @param reverseComplementLeft left read is reverse complement
   * @param templateStartLeft template start of the leftmost read
   * @param reverseComplementRight if the right is reverse complement
   * @param templateStartRight template start of the rightmost read
   */
  public void setValues(int readId, boolean firstRight, boolean reverseComplementLeft, int templateStartLeft, boolean reverseComplementRight, int templateStartRight) {
    mReadId = readId;
    mFirstRight = firstRight;
    mReverseComplementRight = reverseComplementRight;
    mTemplateStartRight = templateStartRight;
    mReverseComplementLeft = reverseComplementLeft;
    mTemplateStartLeft = templateStartLeft;
    mLeftIsVeryPoor = false;
    setAlignments(null, null);
  }

  /**
   * Set data about this hit
   * @param left the alignment result of the leftmost alignment
   * @param right the alignment result of the rightmost alignment
   */
  public void setAlignments(AlignmentResult left, AlignmentResult right) {
    mAlignmentLeft = left;
    mAlignmentScoreLeft = mAlignmentLeft == null ? -1 : left.getScore();
    mAlignmentRight = right;
    mAlignmentScoreRight = mAlignmentRight == null ? -1 : right.getScore();
  }

  public int getAlignmentScoreLeft() {
    return mAlignmentScoreLeft;
  }

  public void setAlignmentScoreLeft(int alignmentScoreLeft) {
    mAlignmentScoreLeft = alignmentScoreLeft;
  }

  public int getAlignmentScoreRight() {
    return mAlignmentScoreRight;
  }

  public void setAlignmentScoreRight(int alignmentScoreRight) {
    mAlignmentScoreRight = alignmentScoreRight;
  }

  public void setLeftVeryPoor(boolean leftIsVeryPoor) {
    mLeftIsVeryPoor = leftIsVeryPoor;
  }

  public int getReadId() {
    return mReadId;
  }

  public AlignmentResult getAlignmentLeft() {
    return mAlignmentLeft;
  }

  public AlignmentResult getAlignmentRight() {
    return mAlignmentRight;
  }

  public boolean isFirstRight() {
    return mFirstRight;
  }

  public int getTemplateStartLeft() {
    return mTemplateStartLeft;
  }

  public boolean isReverseComplementRight() {
    return mReverseComplementRight;
  }

  public int getTemplateStartRight() {
    return mTemplateStartRight;
  }

  public boolean isReverseComplementLeft() {
    return mReverseComplementLeft;
  }
}

